#!/usr/bin/python
# programmer : bbc
# usage:

import sys
import re
import random
import string
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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

def norm(infile,method,flag="single"):
	#For single treatment and single control
  #Check for negative numbers and replace them to zero
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


def main():
	infile = pd.read_table(sys.argv[1],sep="\t")
	#print norm(infile,"upperquartile","single")

if __name__ == '__main__':
	main()
	
