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


def run(infile, designfile, ofile, t, c, collapsemethod):
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
    comparisons = makecomparisions(cfile,treatment_cols,control_cols)
    #Run stattest for each comparision and combine them together
    #df = reformat(cfile, treatment_cols, control_cols, collapsemethod)
    #df.to_csv(ofile,sep="\t",index=False)
    return
  else:
    dfile = pd.read_table(designfile)
    treatment_cols = np.array(t.split(","),dtype=str)
    control_cols = np.array(c.split(","),dtype=str)
    fnames = []
    queuedic = {}
    workers = []
    for tgroup,cgroup in zip(treatment_cols, control_cols):
      #For each comparison group, construct comparison for all replicates
      t_sample = dfile[dfile['group']==tgroup].loc[:,'sample'].unique()
      c_sample = dfile[dfile['group']==cgroup].loc[:,'sample'].unique()
      t_cols = mapcolindex(cfile._get_numeric_data().columns,t_sample)
      c_cols = mapcolindex(cfile._get_numeric_data().columns,c_sample)
      comparisons = makecomparisons(cfile, t_cols, c_cols) 
      cmp_group_name = tgroup+"_vs_"+cgroup
      queuedic[cmp_group_name] = Queue() 
      for comparison_df in comparisons:
        p = Process(target=callstat, args=queuedic[cmp_group_name], comparison_df, outname, multiplier, randomseed)
        workers.append(p)
        p.start()
     for worker in workers:
       worker.join()
      #df = reformat(cfile, t_cols, c_cols,collapsemethod)
      #outname = ofile+"_"+tgroup+"_vs_"+cgroup
      #df.to_csv(outname,sep="\t",index=False)
      #fnames.append(outname)
    return fnames

def callstat(queue, df, outname, multiplier, randomseed):
  queue.put(stattest.runStatinfer(df, outname, multiplier, randomseed))

def makecomparisons(cfile, t_cols, c_cols):
  '''
  Given normalized count file and column indexs of control and treatment,
  generate 4-col dataframe for stattest
  '''
  df_list = []
  for t_index, n_index in zip(t_cols, c_cols):
    df_list.append(cfile.iloc[:,[0,1, t_index+2, c_index+2]])
  if (len(t_cols)>len(c_cols)):#need to treat unpaired treatment samples
    number_of_unpair_treatment = len(t_cols) - len(c_cols)
    averge_control = averageReplicates(cfile, c_cols) 
    for index in range(number_of_unpair_treatment):
      df = cfile.iloc[:,[0,1,index+2]]
      df = df.merge(average_control, on='sgRNA', how='left')
      df = df.fillna(1.0)
      df_list.append(df)
  return df_list

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
  runReformat(args.infile, args.designfile, args.outfile, args.treat, args.ctrl)

def main():
  argparser = prepare_argparser()
  args = argparser.parse_args()
  runReformat(args.infile, args.designfile, args.outfile, args.treat, args.ctrl)
    
if __name__=="__main__":
  main()
