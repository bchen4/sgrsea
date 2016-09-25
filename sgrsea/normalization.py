#!/usr/bin/python
# programmer : bbc
# usage:

import sys
import logging
import numpy as np
import pandas as pd
import argparse as ap
import matplotlib.pyplot as plt
logging.basicConfig(level=10)

def prepare_argparser():
  description = "Normalize count data"
  epilog = "For command line options of each command, type %(prog)% COMMAND -h"
  argparser = ap.ArgumentParser(description=description, epilog = epilog)
  argparser.add_argument("-i","--input",dest = "infile",type=str, required=True, help="input count file matrix")
  argparser.add_argument("--normalize-method",dest="method", default="total", type=str,help ="design file", choices=['total','median','upperquantile'])
  argparser.add_argument("-o","--output",dest = "outfile",type=str,required=True, help="output")
  argparser.add_argument("--split-lib",dest = "splitlib",action='store_true', help="Lib A and B are sequenced separately")
  return argparser

def percentile(n):
  def percentile_(x):
    return np.percentile(x,n)
  percentile_.__name__ = 'percentile_%s' % n
  return percentile_


def norm(infile,method):
  info = infile.loc[:,['sgRNA','Gene']]
  num = infile._get_numeric_data()
  num[num < 0] = 0
  if method == "upperquartile":
    norm_factor = num.apply(percentile(75),axis=0).astype(float)
    smooth_factor = float(norm_factor.mean())  
  elif method == "total":
    norm_factor = num.apply(np.sum,axis=0).astype(float)
    smooth_factor = 10**6 * 1.0
  num_norm = num.div(norm_factor/float(smooth_factor),axis="columns")
  num_norm = num_norm.applymap(lambda x: x+1)
  #num_norm = num_norm.applymap(around)
  norm_df = info.join(num_norm)
  return norm_df

def normalization(infile, outfile, method, splitlib):
  if splitlib: #normalize sub lib separately
    myfile = pd.read_table(infile)
    logging.debug(myfile._get_numeric_data().columns)
    norm_df = pd.DataFrame(columns = ['sgRNA','Gene']+myfile._get_numeric_data().columns.tolist())
    for sublib, group in myfile.groupby('sublib'):
      logging.debug(sublib)
      norm_df = norm_df.append(norm(group, method))
  else:
    norm_df = norm(pd.read_table(infile),method)
  norm_df.to_csv(outfile, sep="\t", index=False)

def main():
  argparser = prepare_argparser()
  args = argparser.parse_args()
  normalization(args.infile, args.outfile, args.method, args.splitlib)
  #infile = pd.read_table(args.infile)
  #logging.debug("read file")
  #print norm(infile,"total")

if __name__ == '__main__':
  main()
  
