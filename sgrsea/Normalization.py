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
  argparser.add_argument("-m","--method",dest="method", default="total", type=str,help ="design file", choices=['total','median','upperquantile'])
  argparser.add_argument("-o","--output",dest = "outfile",type=str,required=True, help="output")
  argparser.add_argument("--split-lib",dest = "splitlib",action='store_true', help="Lib A and B are sequenced separately")



def percentile(n):
	def percentile_(x):
		return np.percentile(x,n)
	percentile_.__name__ = 'percentile_%s' % n
	return percentile_

def around_np(n):# decimal parameter, return integer
	def around_(x):
		return np.around(x,decimals=n)
	around_.__name__= 'around_%s' % n
	return around_

def around(x):
  '''Rounds 0.4 to 0 and 0.5 to 1. Only deals with positive numbers'''
  try:
    x = float(x)
  except:
    logging.error("There are values "+str(x)+" cannot be converted to float. Please make sure your count data start from the 3rd column of your data frame.")
    sys.exit(1)
  else:
    if x - int(x) <= 0.4:
      return int(x)
    else:
      return int(x)+1

def norm(infile,method):
  num = infile._get_numeric_data()
  num[num < 0] = 0
  if method == "upperquartile":
		norm_factor = infile.iloc[:,2:].apply(percentile(75),axis=0).astype(float)
		smooth_factor = float(norm_factor.mean())	
  elif method == "total":
		norm_factor = infile.iloc[:,2:].apply(np.sum,axis=0).astype(float)
		smooth_factor = 10**6 * 1.0
	#Normalize for sample size	
  infile.iloc[:,2:] = infile.iloc[:,2:].div(norm_factor/float(smooth_factor),axis="columns")
  infile.iloc[:,2:] = infile.iloc[:,2:].applymap(lambda x: x+1)

  infile.iloc[:,2:] = infile.iloc[:,2:].applymap(around)
  return infile

def normalization(infile, outfile, method, splitlib):
  if splitlib: #normalize sub lib separately
    norm_df = pd.DataFrame(columns = infile.columns)
    for sublib, group in infile.groupby('sublib'):
      norm_df.append(norm(infile, method))
  else:
    norm_df = norm(infile,method)
  norm_df.to_csv(outfile, sep="\t", index=False)

def main():
	argparser = prepare_argparser()
  args = argparser.parse_args()
  normalization(args.infile, args.outfile, args.method, args.splitlib)
  #infile = pd.read_table(sys.argv[1],sep="\t")
	#print norm(infile,"upperquartile","single")

if __name__ == '__main__':
	main()
	
