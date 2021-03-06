# programmer : bbc
# usage:

import sys
import argparse as ap
import logging
import pandas as pd
import numpy as np
import stattest
import reformatCountTable
from multiprocessing import Pool
from scipy.stats.mstats import gmean

logging.basicConfig(level=10)


def prepare_argparser():
  description = "Rank gene"
  epilog = "For command line options of each command, type %(prog)% COMMAND -h"
  argparser = ap.ArgumentParser(description=description, epilog = epilog)
  argparser.add_argument("-i","--input",dest = "infile",type=str,required=True, help="input BAM file")
  argparser.add_argument("-o","--output",dest = "outfile",type=str,required=True, help="output")
  argparser.add_argument("-d","--design",dest = "designfile",type=str, help="design file")
  argparser.add_argument("-t","--treatment",dest="treat",type=str,required=True, help = "columns/name of treatment samples")
  argparser.add_argument("-c","--control",dest="ctrl",type=str,required=True, help="columns/name of control samples")
  argparser.add_argument("--multiplier",dest="multiplier",type=int,default=50, help="Multiplier to generate background")
  argparser.add_argument("--random-seed",dest="randomSeed",type=int,default=None, help="Random seed to control permutation process")
  argparser.add_argument("--collapse-replicates", dest="collapsemethod", type=str, help = "Way to collapse replicates", default="None", choices=['auto','stack','mean','None'])
  return(argparser)


def generank(infile, designfile, ofile, t, c, collapse="None", multiplier=50, randomseed=0):
  if collapse == "None":
    gmrank(infile, designfile, ofile, t, c, multiplier, randomseed)
  else:
    files = reformatCountTable.runReformat(infile, designfile, ofile, t,c, collapse)
    if len(files)==0:
      logging.error("There is no input files for sgRSEA to run. Exit")
      sys.exit(1)
    else:
      rungroupcmp(files, ofile, multiplier, randomseed, 5, False)


def gmrank(infile, designfile, ofile, t, c, multiplier=50, randomseed=0):
  '''
  First construct comparisons:
   1. If treatment and control are numeric, there is no need of design file.
      Each column will be treated as replicates;
   2. If treatment and control are strings, they are assumed to be group. Each
      pair of replicates in that group will be treated separately;
  For treatment replicates don't have matching control, the average of all 
  controls will be calculated and comparisons will be made.

  '''
  try:
    cfile = pd.read_table(infile)
  except IOError,message:
    print >> sys.stderr, "cannot open file",message
    sys.exit(1)
  #logging.debug(t)
  #logging.debug(c)
  try:
    treatment_cols = np.array(t.split(","), dtype=int)
    control_cols = np.array(c.split(","),dtype=int)
    #comparisons = makecomparisons(cfile,treatment_cols,control_cols,ofile)
    #Run stattest for each comparision and combine them together
    #rungroupcmp(comparisons, ofile, multiplier, randomseed)
  except:# Use design file
    dfile = pd.read_table(designfile)
    treatment_cols = np.array(t.split(","),dtype=str)
    control_cols = np.array(c.split(","),dtype=str)
    fnames = []
    for tgroup,cgroup in zip(treatment_cols, control_cols):
      #For each comparison group, construct comparisons for all replicates
      #logging.debug(tgroup)
      #logging.debug(cgroup)
      t_sample = dfile[dfile['group']==tgroup].loc[:,'sample'].unique()
      c_sample = dfile[dfile['group']==cgroup].loc[:,'sample'].unique()
      t_cols = mapcolindex(cfile._get_numeric_data().columns,t_sample)
      c_cols = mapcolindex(cfile._get_numeric_data().columns,c_sample)
      logging.debug(t_cols)
      logging.debug(c_cols)
      comparisons = makecomparisons(cfile, t_cols, c_cols, ofile) 
      rungroupcmp(comparisons, ofile+"_"+tgroup+"_vs_"+cgroup, multiplier, randomseed)  
  else:
    comparisons = makecomparisons(cfile,treatment_cols,control_cols,ofile)
    #Run stattest for each comparision and combine them together
    rungroupcmp(comparisons, ofile, multiplier, randomseed)



def rungroupcmp(inputdfnamelist, outprefix, multiplier, randomseed, number_of_workers=5,gm_flag=True):
  '''
  Run all comparisons for a group
  '''
  work_num = min(len(inputdfnamelist), number_of_workers)
  work_pool = Pool(work_num)
  #make output name list
  outnamelist = []
  for name in inputdfnamelist:
    #logging.debug(name)
    outnamelist.append(name.replace(".forStat",""))
  arglist = zip(inputdfnamelist,outnamelist,[multiplier]*len(inputdfnamelist), [randomseed]*len(inputdfnamelist) )
  resultList = work_pool.map(callstat, arglist)
  work_pool.close()
  work_pool.join()
#get result sgRSEA result file names
  ##logging.debug(resultList)
  if len(resultList)>1 and gm_flag:#merge results together and calculate geometric mean
  #change df column names and merge
    pos_rank_cols = []#record new rank column names for future use
    neg_rank_cols = []
    ini_dfname = resultList[0]
    res_df = pd.read_table(ini_dfname)
    res_df = res_df.iloc[:,[0,1,7,8]]
    #logging.debug(res_df.columns)
    pos_rank_cols.append(res_df.columns.tolist()[2])
    neg_rank_cols.append(res_df.columns.tolist()[3])
    #logging.debug(pos_rank_cols)
    #logging.debug(neg_rank_cols)
    for dfname in resultList[1:]:
      df = pd.read_table(dfname)
      df = df.iloc[:,[0,7,8]]
      pos_rank_cols.append(df.columns.tolist()[1])
      neg_rank_cols.append(df.columns.tolist()[2])
      
      res_df = res_df.merge(df, on=['Gene'])
  #Calculate geometric mean for the ranks
  #logging.debug("pos_rank_cols and neg_rank_cols")
  #logging.debug(pos_rank_cols)
  #logging.debug(neg_rank_cols)

    res_df['pos_geomean'] = gmean(res_df.loc[:, pos_rank_cols], axis=1)
    res_df['neg_geomean'] = gmean(res_df.loc[:, neg_rank_cols], axis=1)
    res_df = res_df.sort_values(by='pos_geomean', ascending=True)
    res_df = res_df.reset_index(drop=True)
    res_df['pos_geomean_rank'] = res_df.index + 1
    res_df['neg_geomean_rank'] = res_df.shape[0] - res_df['pos_geomean_rank']

    res_df.to_csv(outprefix+".gm.sgRSEA.xls", sep="\t", index=False)


def callstat(args):
  return stattest.runStatinfer(*args)

def makecomparisons(cfile, t_cols, c_cols, outprefix):
  '''
  Given normalized count file and column indexs of control and treatment,
  generate 4-col dataframe for stattest, save to file and return file name
  list
  '''
  df_name_list = []
  for t_index, c_index in zip(t_cols, c_cols):
    new_df = cfile.iloc[:,[0,1, t_index+2, c_index+2]]
    new_df_name = "_".join([outprefix, cfile.columns.tolist()[t_index+2],"vs", cfile.columns.tolist()[c_index+2]])+".forStat"
    df_name_list.append(new_df_name)
    new_df.to_csv(new_df_name,sep="\t",index=False)
  if (len(t_cols)>len(c_cols)):#need to treat unpaired treatment samples
    average_control = averageReplicates(cfile, c_cols) 
    for index in range(len(c_cols),len(t_cols)):
      df = cfile.iloc[:,[0,1,index+2]]
      df = df.merge(average_control, on='sgRNA', how='left')
      df = df.fillna(1.0)
      df_name = "_".join([outprefix, cfile.columns.tolist()[index+2],"vs","averageCtrl"])+".forStat"
      df_name_list.append(df_name)
      df.to_csv(df_name,sep="\t",index=False)
  return df_name_list

def mapcolindex(header,samples):
  col = []
  for i,v in enumerate(header):
    if v in (samples):
      col.append(i)
  return col

def averageReplicates(cfile, cols):
  nfile = cfile._get_numeric_data()
  cfile['ave_ctrl'] = nfile.iloc[:,cols].mean(axis=1)
  return cfile.loc[:,['sgRNA','ave_ctrl']]

def run(args):
  generank(args.infile, args.designfile, args.outfile, args.treat, args.ctrl, args.collapsemethod, args.multiplier, args.randomSeed)


def main():
  argparser = prepare_argparser()
  args = argparser.parse_args()
  run(args.infile, args.designfile, args.outfile, args.treat, args.ctrl, args.collapsemethod, args.multiplier, args.randomSeed)
    
if __name__=="__main__":
  main()
