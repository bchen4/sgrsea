#!/usr/bin/python
# programmer : bbc
# usage:

import sys
import random
import logging
import copy
import math
import numpy as np
import pandas as pd
import argparse as ap
from statsmodels.stats.multitest import multipletests

logging.basicConfig(format='%(levelname)s:%(message)s',level=logging.DEBUG)

def prepare_argparser():
  description = "sgRSEA main stat function"
  epilog "For command line options of each command, type %(prog)% COMMAND -h"
  argparser = ap.ArgumentParser(description=description, epilog = epilog)
  argparser.add_argument("-i","--input",dest = "infile",type=str,required=True, help = "sgRSEA input file, 4 columns")
  argparser.add_argument("-o","--output",dest = "outfile",type=str,required=True, help = "output file name")
  argparser.add_argument("-m","--multiplier",dest = "multiplier",type=int, default = 30,required=True, help = "Multiplier to generate background")
  argparser.add_argument("--bgtag",dest = "bgtag",type=str, help = "Sting to identify control sgRNAs")
  argparser.add_argument("--bg-row-start",dest = "bgrowstart",type=int,default = -1, help = "Row count of the start of control sgRNA block")
  argparser.add_argument("--bg-row-stop",dest = "bgrowstart",type=int, default=-1, help = "Row count of the stop of control sgRNA block")
  return(argparser)

def getBackground(infile,nontag="",tagStart=0,tagStop=0):
  '''Get background data frame from either string tag or row range.
  '''
  if len(nontag)>0:
    bgFile = infile[infile.iloc[:,0].str.contains(nontag)]
    dataFile = infile[infile.iloc[:,0].str.contains(nontag)==False]
  elif tagStart<tagStop and tagStart<infile.shape[0]: #use the row range to get file
    tagStart -= 1
    tagStop = min(infile.shape[0],tagStop)
    bgFile = infile.iloc[tagStart:tagStop,:]
    dataFile = infile.iloc[0:tagStart-1,:].append(infile.iloc[tagStop:,:])
  else:#Non background
    logging.info("No designed NonTargeting sgRNAs, use input data as backgound.")
    bgFile = copy.copy(infile)
    dataFile = infile

  logging.debug("bgFile: "+str(bgFile.shape[0]))
  logging.debug("dataFile: "+str(dataFile.shape[0]))
  if bgFile.shape[0]==0: #backgound file is empty
    logging.info("There is no NonTargeting sgRNA found. Use input data as background.")
    bgFile = copy.copy(infile)
    dataFile = infile
  return (dataFile,bgFile)

def addZstat(data_df, pNull):
  df = copy.copy(data_df)
  df['dnorm'] =  [math.sqrt(pNull * (1.0 - pNull))] * df.shape[0]
  df['total'] = (df.loc[:,'treat_1'] + df.loc[:,'ctrl_1']).values
  df['t_sqrt'] = df.loc[:,'total'].apply(lambda x: math.sqrt(x))
  df['pmme'] = (df.loc[:,'treat_1']/df.loc[:,'total']).values
  df['pnull'] = [pNull] * df.shape[0]
  df['zstat'] = (df.loc[:,'pmme'] - df.loc[:,'pnull']) * df.loc[:,'t_sqrt'] / df.loc[:,'dnorm']
  df['maxmean'] = [0]*df.shape[0]
  df['sgcount'] = [0]*df.shape[0]
  return df.loc[:,['sgRNA','gene','treat_1','ctrl_1','zstat','maxmean','sgcount']]

def dataFilter(dataFile,sg_min = 1):
  '''Get table of sgRNA number per gene distribution. Filter out the genes with sgRNA number less than sg_min'''
  filter_data = dataFile.groupby('gene').filter(lambda x : len(x)>=sg_min)
  logging.debug(type(filter_data))
  logging.debug(filter_data.columns.values)
  #filter_data['zstat'] = filter_data.loc[:,['treat_1','ctrl_1']].apply(lambda row: zStat(row.to_frame(['treat_1','ctrl_1']),0.4))
  gene_label = filter_data.loc[:,'gene']
  group_data = filter_data.groupby('gene')
  return (group_data,gene_label)

def getNewHeader(treatments,controls):
  if len(np.intersect1d(treatments,controls))>0:
    return None
  header = ['sgRNA','gene']
  for i in range(len(treatments)+len(controls)):
    header.append('')
  count = 0
  for i in treatments:
    header[int(i)-1] = 'treat_'+str(count+1)
    count += 1
  count = 0
  for j in controls:
    header[int(j)-1] = 'ctrl_'+str(count+1)
    count += 1
  return header

def pMME(treat,ctrl):
  '''Input are two pd.Series'''
  sum_treat = treat.sum().item()
  sum_ctrl = ctrl.sum().item()
  sum_all = sum_treat + sum_ctrl
  try:
    pmme = sum_treat * 1.0 / sum_all
  except ZeroDivisionError:
    #logging.warning("ZeroDivisionError in pMME calculation, which indicates no reads count for treatment. Return 0 instead. Continuing")
    pmme = 0.0
  if not (0.0<=pmme<=1.0):
    logging.error("pMME range error. Exit")
    sys.exit(1)
  return pmme

def zStat(sg_row,pNull):
  if pNull == 0 or pNull ==1:
    return None
  else:
    zscore = (pMME(sg_row.treat_1,sg_row.ctrl_1)-pNull) * math.sqrt(sg_row.sum(0).sum()) / math.sqrt(pNull*(1-pNull)) 
    return zscore

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
  logging.debug(maxmean_df.head(10))
  return maxmean_df

def sampleNull(genelist,bgFile):
  sample_size = len(genelist)
  samples = bgFile.iloc[np.random.randint(0,len(bgFile),size=sample_size)]
  #make samples into gene_matrix structure list
  samples['gene'] = genelist.values
  logging.debug(samples.head(10))
  group_sample = samples.groupby('gene')
  return group_sample


def maxmeanSampleNull(treat_stat,bgFile,pNull,multiplier=10):
  '''Call sampleNull() multiplier times. Only keep maxmean dataframe to save memory.'''
  null_maxmean = pd.DataFrame()
  np.random.seed(8512)
  for i in range(multiplier):
    null_group = sampleNull(treat_stat,bgFile)
    newdf = getMatrixMaxmean(null_group,pNull)
    null_maxmean = null_maxmean.append(newdf)
  logging.debug("Null_maxmean_df")
  logging.debug(null_maxmean[null_maxmean['sgcount'].isin([8,9])])
  return null_maxmean

def standardizeFactor(null_maxmean):
  '''Input: null_maxmean_df, Returns a df with 'sgcount', 'mean','std' '''
  group_null_maxmean = null_maxmean.groupby('sgcount')
  factor_df = pd.DataFrame({'mean':group_null_maxmean.maxmean.agg(np.mean),'std':group_null_maxmean.maxmean.agg(np.std)}).reset_index()
  return factor_df
  
def standardize(datarow,factors):
  '''Datarow is a pd.Series, index is ['gene','maxmean','sgcount']'''
  factor = factors[factors['sgcount']==datarow.sgcount]
  if factor.shape[0]==1:#Should be only one record
    f_mean = factor.iloc[0]['mean']
    f_std = factor.iloc[0]['std']
  elif factor.shape[0]==0:
    logging.warning("There is no record for genes with "+str(datarow.sgcount+" sgRNAs. "))
  else:#more than 1 record
    logging.warning("There are more than 1 record for sgRNA count "+str(datarow.sgcount+" group. Check data."))
  smm = (datarow.maxmean-f_mean)*1.0/f_std
  return smm

def standardizeDF(maxmean_df,sFactors):
  maxmean_df['sMaxmean'] = maxmean_df.apply(lambda x: standardize(x,sFactors),axis=1)
  return maxmean_df
  
def pvalue(tn,smm_null):
  '''smm_null is a dataframe'''
  p = (smm_null[smm_null['sMaxmean']>tn].shape[0]+1.0)/smm_null.shape[0]
  return p

def getPQ(data_maxmean_std,null_maxmean_std):
  data_maxmean_std['pval'] = data_maxmean_std.apply(lambda x: pvalue(x['sMaxmean'],null_maxmean_std),axis=1)
  data_maxmean_std['qval'] = multipletests(data_maxmean_std.loc[:,'pval'],method='fdr_bh')[1]
  return data_maxmean_std

def runStatinfer(infile,nontag,tagStart,tagStop,multiplier):
  logging.info("Start to run test.")
  #reset dataframe column names
  old_header = list(infile.columns.values)
  new_header = getNewHeader(treatments,controls)
  infile.columns = new_header
  (dataFile,bgFile) = getBackground(infile,nontag,tagStart,tagStop)
  p0 = pMME(bgFile.loc[:,['treat_1']],bgFile.loc[:,['ctrl_1']])
  if p0 ==0 or p0 ==1:
    logging.error("pMME for background equals to 0/1, indicating no counts for treatment or control. Please check your data. Exit.")
    sys.exit(1)
  dataFile = addZstat(dataFile,p0)
  #logging.debug(dataFile.head(10))
  #dataFile.to_csv("test_real_zscore.txt",sep="\t",header=True,index=False)
  bgFile = addZstat(bgFile, p0)
  (treat_group,genelist) = dataFilter(dataFile,1)
  #dataStat = treat_group.size()
  #logging.info("Calculating treatment maxmean...")
  data_maxmean_df = getMatrixMaxmean(treat_group)
  #data_maxmean_df.to_csv("test_real_maxmean.txt",sep="\t",header=True,index=False)
  #logging.info("Sampling null distribution...")
  null_maxmean_df = maxmeanSampleNull(genelist,bgFile,p0,multiplier)
  #logging.debug("Get standardize factors")
  factor_df = standardizeFactor(data_maxmean_df)
  logging.debug(factor_df)
  #logging.info("Standardization...")
  data_sdf = standardizeDF(data_maxmean_df,factor_df)
  null_sdf = standardizedDF(null_maxmean_df,factor_df)
  #sDf.to_csv("test_real_standardize.txt",sep="\t",header=True,index=False)
  fdf = getPQ(data_sdf,null_sdf)
  return fdf

def main():
  argparser = prepare_argparser()
  args = argparser.parse_args()
  runStatinfer(args.infile,args.outfile,args.bgtag, args.bgstart, args.bgstop, args.multiplier)
  
if __name__ == '__main__':
	main()
	
