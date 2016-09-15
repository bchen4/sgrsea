#!/usr/bin/python
# programmer : bbc
# usage:

import sys
import argparse as ap
import logging
import pandas as pd
import numpy as np

logging.basicConfig(level=10)


def prepare_argparser():
  description = "Reformat the count table to sgRSEA format"
  epilog = "For command line options of each command, type %(prog)% COMMAND -h"
  argparser = ap.ArgumentParser(description=description, epilog = epilog)
  argparser.add_argument("-i","--input",dest = "infile",type=str,required=True, help="input BAM file")
  argparser.add_argument("-o","--output",dest = "outfile",type=str,required=True, help="output")
  argparser.add_argument("-d","--design",dest = "designfile",type=str, help="output")
  argparser.add_argument("-t","--treatment",dest="treat",type=str,required=True, help = "columns/name of treatment samples")
  argparser.add_argument("-c","--control",dest="ctrl",type=str,required=True, help="columns/name of control samples")
  argparser.add_argument("--t-lable",dest="treatlabel",type=str, help = "label of treatment samples")
  argparser.add_argument("--c-label",dest="ctrllabel",type=str, help="label of control samples")
  return(argparser)


def runReformat(infile, designfile, ofile, t, c):
  '''
  If treatment and control are numeric, there is no need of design file.
  '''
  try:
    cfile = pd.read_table(infile)
  except IOError,message:
    print >> sys.stderr, "cannot open file",message
    sys.exit(1)
  if not designfile:
    treatment_cols = np.array(t.split(","),dtype=int)
    control_cols = np.array(c.split(","),dtype=int)
    df = reformat(cfile, treatment_cols, control_cols)
    df.to_csv(ofile,sep="\t",index=False)
    return [ofile]
  else:
    dfile = pd.read_table(designfile)
    treatment_cols = np.array(t.split(","),dtype=str)
    control_cols = np.array(c.split(","),dtype=str)
    fnames = []
    for tgroup,cgroup in zip(treatment_cols, control_cols):
      t_label = dfile[dfile['group']==tgroup].loc[:,'label'].unique()
      c_label = dfile[dfile['group']==cgroup].loc[:,'label'].unique()
      t_cols = mapcolindex(cfile._get_numeric_data().columns,t_label)
      c_cols = mapcolindex(cfile._get_numeric_data().columns,c_label)
      df = reformat(cfile, t_cols, c_cols)
      outname = ofile+"_"+t_label[0]+"_vs_"+c_label[0]
      df.to_csv(outname,sep="\t",index=False)
      fnames.append(outname)
    return fnames


def mapcolindex(header,labels):
  col = []
  for i,v in enumerate(header):
    if v in (labels):
      col.append(i)
  return col

def reformat(cfile, treatcols, ctrlcols):
  nfile = cfile._get_numeric_data()
  logging.debug(nfile.shape)
  if len(treatcols)!=len(ctrlcols):#take average
    cfile['treat_ave'] = nfile.iloc[:,treatcols].mean(axis=1)
    cfile['ctrl_ave'] = nfile.iloc[:,ctrlcols].mean(axis=1)
    df = cfile.loc[:,['sgRNA','Gene','treat_ave','ctrl_ave']]
  else:#combine
    df = cfile.loc[:,['sgRNA','Gene']].join(nfile.iloc[:,[treatcols[0],ctrlcols[0]]])
    df.columns = ['sgRNA','Gene','treatment','control']
    for tc, cc in zip(treatcols[1:],ctrlcols[1:]):
      logging.debug(tc)
      logging.debug(cc)
      newdf = cfile.loc[:,['sgRNA','Gene']].join(nfile.iloc[:,[tc,cc]])
      newdf.columns = df.columns
      df = df.append(newdf)
  return df


def main():
  argparser = prepare_argparser()
  args = argparser.parse_args()
  runReformat(args.infile, args.designfile, args.outfile, args.treat, args.ctrl)
    
if __name__=="__main__":
  main()