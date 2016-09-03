#!/usr/bin/python
# programmer : bbc
# usage:

import sys
import re
import random
import string
import logging
import copy
import math
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from datetime import datetime
from multiprocessing import Process, Queue

logging.basicConfig(format='%(levelname)s:%(message)s',level=logging.DEBUG)

def maxMean(gene_df):
  scores = gene_df.loc[:,'zstat']
  sc_po = sum(scores[scores>=0])*1.0
  sc_ne = sum(scores[scores<0])*1.0
  if abs(sc_po) >= abs(sc_ne):
    return sc_po/gene_df.shape[0]
  else:
    return sc_ne/gene_df.shape[0]

def getMatrixMaxmean(gene_group):
  '''Gene group is a groupby object.
     This function returns a dataframe with 'gene', 'maxmean', 'sgCount'    
  '''
  genes = []
  maxmeans = []
  sgCount = []
  for g, gdf in gene_group:
    maxmeans.append(maxMean(gdf))
    genes.append(g)
    sgCount.append(gdf.shape[0])
  #logging.debug(maxmean_df.head(10))
      #sgCount.append(gene_matrix[g].shape[0])
  maxmean_df = pd.DataFrame({'gene':genes,'maxmean':maxmeans,'sgcount':sgCount})
  #logging.debug(maxmean_df.head(10))
  return maxmean_df

def worker(infile_group, out_q):
  try:
    outdf = getMatrixMaxmean(infile_group)
  except:
    out_q.put(pd.DataFrame())
  else:
    out_q.put(outdf)

def multi(infile_group,multiplier):
  
  out_q = Queue()
  procs = []
  for i in range(multiplier):
    p = Process(target=worker,args=(infile_group,out_q))
    procs.append(p)
    p.start()

  resultdict = []
  logging.debug("Job submitted")
  #out_q.close()
  #out_q.join()
  for p in procs:
    logging.debug(p.pid)
    p.join()
  logging.debug("processed finished")
  for i in range(multiplier):
    resultdict.append(out_q.get())

  return resultdict

def sample_multi(infile_group,multiplier):
  result = multi(infile_group,multiplier)
  df = pd.DataFrame()
  for v in result:
    print type(v)
    df = df.append(v)
  print (df.shape)
  print (df.head(10))

def testGroup(infile):
  infile_group = infile.groupby('gene')
  start = datetime.now()
  #getMatrixMaxmean(infile_group)
  #stop = datetime.now()
  #logging.debug(stop-start)
  #start = stop
  sample_multi(infile_group,3)
  stop = datetime.now()
  logging.debug(stop - start)

def main():
  infile = pd.read_table(sys.argv[1],sep="\t")
  infile = infile.iloc[0:1000,:]
  mdf = testGroup(infile)
  
if __name__ == '__main__':
	main()
	
