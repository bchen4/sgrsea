# programmer : bbc
# usage:

import sys
import argparse as ap
import logging
import pandas as pd
import numpy as np

logging.basicConfig(level=10)


def prepare_argparser():
  description = "Rank gene"
  epilog = "For command line options of each command, type %(prog)% COMMAND -h"
  argparser = ap.ArgumentParser(description=description, epilog = epilog)
  argparser.add_argument("-i","--input",dest = "infile",type=str,required=True, help="input BAM file")
  argparser.add_argument("-o","--output",dest = "outfile",type=str,required=True, help="output")
  argparser.add_argument("-d","--design",dest = "designfile",type=str, help="output")
  argparser.add_argument("-t","--treatment",dest="treat",type=str,required=True, help = "columns/name of treatment samples")
  argparser.add_argument("-c","--control",dest="ctrl",type=str,required=True, help="columns/name of control samples")
  return(argparser)


def run(infile, designfile, ofile, t, c, multiplier=50, randomseed=0):
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
  if not designfile:
    treatment_cols = np.array(t.split(","),dtype=int)
    control_cols = np.array(c.split(","),dtype=int)
    comparisons = makecomparisions(cfile,treatment_cols,control_cols,outprefix)
    #Run stattest for each comparision and combine them together
    rungroupcmp(comparisons, ofile, multiplier, randomseed)
    #df = reformat(cfile, treatment_cols, control_cols, collapsemethod)
    #df.to_csv(ofile,sep="\t",index=False)
  else:
    dfile = pd.read_table(designfile)
    treatment_cols = np.array(t.split(","),dtype=str)
    control_cols = np.array(c.split(","),dtype=str)
    fnames = []
    queuedic = {}
    workers = []
    for tgroup,cgroup in zip(treatment_cols, control_cols):
      #For each comparison group, construct comparisons for all replicates
      t_sample = dfile[dfile['group']==tgroup].loc[:,'sample'].unique()
      c_sample = dfile[dfile['group']==cgroup].loc[:,'sample'].unique()
      t_cols = mapcolindex(cfile._get_numeric_data().columns,t_sample)
      c_cols = mapcolindex(cfile._get_numeric_data().columns,c_sample)
      comparisons = makecomparisons(cfile, t_cols, c_cols, outprefix) 
      
      
def rungroupcmp(inputdflist, outname, multiplier, randomseed, number_of_workers=5):
  '''
  Run all comparisons for a group
  '''
  work_num = min(len(inputdflist), number_of_workers)
  work_pool = Pool(work_num)
  for comparison_df in comparisons:
    p = Process(target=callstat, args=queuedic[cmp_group_name], comparison_df, outname, multiplier, randomseed)
    workers.append(p)
    p.start()
  for worker in workers:
    worker.join()

def callstat(queue, df, outname, multiplier, randomseed):
  queue.put(stattest.runStatinfer(df, outname, multiplier, randomseed))

def makecomparisons(cfile, t_cols, c_cols, outprefix):
  '''
  Given normalized count file and column indexs of control and treatment,
  generate 4-col dataframe for stattest, save to file and return file name
  list
  '''
  df_name_list = []
  for t_index, n_index in zip(t_cols, c_cols):
    new_df = cfile.iloc[:,[0,1, t_index+2, c_index+2]]
    new_df_name = "_".join([outprefix, str(t_index),str(c_index)])+".forStat.txt"
    df_name_list.append(new_df_name)
    new_df.to_csv(new_df_name,sep="\t",index=False)
  if (len(t_cols)>len(c_cols)):#need to treat unpaired treatment samples
    number_of_unpair_treatment = len(t_cols) - len(c_cols)
    averge_control = averageReplicates(cfile, c_cols) 
    for index in range(number_of_unpair_treatment):
      df = cfile.iloc[:,[0,1,index+2]]
      df = df.merge(average_control, on='sgRNA', how='left')
      df = df.fillna(1.0)
      df_name = "_".join([outprefix, str(c_cols+index),"averageCtrl"])+".forStat.txt"
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



def main():
  argparser = prepare_argparser()
  args = argparser.parse_args()
  run(args.infile, args.designfile, args.outfile, args.treat, args.ctrl, args.multiplier, args.randomseed)
    
if __name__=="__main__":
  main()
