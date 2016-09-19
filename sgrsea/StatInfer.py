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
#from statsmodels.stats.multitest import multipletests

logging.basicConfig(format='%(levelname)s:%(message)s',level=logging.DEBUG)

def prepare_argparser():
  description = "sgRSEA main stat function"
  epilog = "For command line options of each command, type %(prog)% COMMAND -h"
  argparser = ap.ArgumentParser(description=description, epilog = epilog)
  argparser.add_argument("-i","--input",dest = "infile",type=str,required=True, help = "sgRSEA input file, 4 columns")
  argparser.add_argument("-o","--output",dest = "outfile",type=str,required=True, help = "output file name")
  argparser.add_argument("-m","--multiplier",dest = "multiplier",type=int, default = 30,required=True, help = "Multiplier to generate background")
  argparser.add_argument("--bgtag",dest = "bgtag",type=str, default = "",help = "Sting to identify control sgRNAs")
  argparser.add_argument("--bg-row-start",dest = "bgrowstart",type=int,default = -1, help = "Row count of the start of control sgRNA block")
  argparser.add_argument("--bg-row-stop",dest = "bgrowstop",type=int, default=-1, help = "Row count of the stop of control sgRNA block")
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
  df['zstat'] = df.apply(lambda x: zStat(x,pNull),axis=0)
  return df

def dataFilter(dataFile,sg_min = 1):
  '''Get table of sgRNA number per gene distribution. Filter out the genes with sgRNA number less than sg_min'''
  filter_data = dataFile.groupby('Gene').filter(lambda x : len(x)>=sg_min)
  logging.debug(type(filter_data))
  logging.debug(filter_data.columns.values)
  #filter_data['zstat'] = filter_data.loc[:,['treat_1','ctrl_1']].apply(lambda row: zStat(row.to_frame(['treat_1','ctrl_1']),0.4))
  gene_label = filter_data.loc[:,'Gene']
  return (filter_data,gene_label)


def pMME(df):
  '''1st numeric col is treatment and 2nd numeric col is control'''
  ndf = df._get_numeric_data()
  sum_col = ndf.sum(axis=1)
  try:
    pmme = sum_col[0].item() * 1.0 / sum_col.sum()
  except ZeroDivisionError:
    #logging.warning("ZeroDivisionError in pMME calculation, which indicates no reads count for treatment. Return 0 instead. Continuing")
    pmme = -1
  if not (0.0<=pmme<=1.0):
    logging.error("pMME range error. Exit")
    sys.exit(1)
  return pmme

def zStat(sg_row,pNull):
  if pNull == 0 or pNull ==1:
    return None
  else:
    zscore = (pMME(sg_row)-pNull) * math.sqrt(sg_row.sum(0).sum()) / math.sqrt(pNull*(1-pNull)) 
    return zscore

def maxMean(gene_zstat):
  scores = gene_zstat
  sc_po = sum(scores[scores>=0])*1.0
  sc_ne = sum(scores[scores<0])*1.0
  if abs(sc_po) >= abs(sc_ne):
    return sc_po/len(scores)
  else:
    return sc_ne/len(scores)

def getMatrixMaxmean(filtered_zdf):
  '''
     This function returns a dataframe with 'Gene', 'maxmean', 'sgCount'    
  '''
  maxmean_df = filtered_zdf.groupby("Gene").agg({'zstat':maxMean,'sgRNA':'count'}).reset_index()
  maxmean_df.columns = ["Gene","maxmean","sgCount"]
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


def maxmeanSampleNull(genelist,bgFile,multiplier=10):
  '''Call sampleNull() multiplier times. Only keep maxmean dataframe to save memory.'''
  null_maxmean = pd.DataFrame()
  np.random.seed(8512)
  for i in range(multiplier):
    null_group = sampleNull(genelist,bgFile)
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
  pos_p = (smm_null[smm_null['sMaxmean']>tn].shape[0]+1.0)/(1.0+smm_null.shape[0])
  return pos_p

def getPQ(data_maxmean_std,null_maxmean_std):
  data_maxmean_std['pos_p'] = data_maxmean_std.apply(lambda x: pvalue(x['sMaxmean'],null_maxmean_std),axis=1)
  data_maxmean_std['neg_p'] = 1.0 -data_maxmean_std['pos_p'] 
  data_maxmean_std['pos_q'] = multipletests(data_maxmean_std.loc[:,'pos_p'],method='fdr_bh')[1]
  data_maxmean_std['neg_q'] = multipletests(data_maxmean_std.loc[:,'neg_p'],method='fdr_bh')[1]
  data_maxmean_std['pos_rank'] = data_maxmean_std['pos_p'].rank()
  data_maxmean_std['neg_rank'] = data_maxmean_std['neg_p'].rank()
  return data_maxmean_std

def runStatinfer(infile,outfile,nontag,tagStart,tagStop,multiplier):
  #reset dataframe column names
  #old_header = list(infile.columns.values)
  #new_header = getNewHeader(treatments,controls)
  #infile.columns = new_header
  if (len(nontag)!=0) or (tagStart>0 and tagStop>0):
    (dataFile,bgFile) = getBackground(infile,nontag,tagStart,tagStop)
    p0 = pMME(bgFile)
  else:#Use dataset as background
    p0 = pMME(dataFile)
  if p0 ==0 or p0 ==1:
    logging.error("pMME for background equals to 0/1, indicating no counts for treatment or control. Please check your data. Exit.")
    sys.exit(1)
  
  dataFile = addZstat(dataFile,p0)
  #logging.debug(dataFile.head(10))
  #dataFile.to_csv("test_real_zscore.txt",sep="\t",header=True,index=False)
  (filtered_data,genelist) = dataFilter(dataFile,1)
  #dataStat = treat_group.size()
  #logging.info("Calculating treatment maxmean...")
  data_maxmean_df = getMatrixMaxmean(filtered_data)
  #data_maxmean_df.to_csv("test_real_maxmean.txt",sep="\t",header=True,index=False)
  #logging.info("Sampling null distribution...")
  if isinstance(bgFile, pd.DataFrame): 
    bgFile = addZstat(bgFile, p0)
  else:
    bgFile = filtered_data
  null_maxmean_df = maxmeanSampleNull(genelist,bgFile,multiplier)
  #logging.debug("Get standardize factors")
  factor_df = standardizeFactor(null_maxmean_df)
  logging.debug(factor_df)
  #logging.info("Standardization...")
  data_sdf = standardizeDF(data_maxmean_df,factor_df)
  null_sdf = standardizedDF(null_maxmean_df,factor_df)
  #sDf.to_csv("test_real_standardize.txt",sep="\t",header=True,index=False)
  fdf = getPQ(data_sdf,null_sdf)
  fdf.to_csv(outfile,sep="\t",index=False)

def main():
  argparser = prepare_argparser()
  args = argparser.parse_args()
  infile = pd.read_table(args.infile)
  #df = getMatrixMaxmean(infile)
  #df.to_csv(args.outfile,sep="\t",index=False)
  runStatinfer(args.infile,args.outfile,args.bgtag, args.bgrowstart, args.bgrowstop, args.multiplier)
  
if __name__ == '__main__':
	main()
	
