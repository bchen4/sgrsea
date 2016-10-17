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
import datetime
from statsmodels.stats.multitest import multipletests

logging.basicConfig(format='%(levelname)s:%(message)s',level=logging.DEBUG)

def prepare_argparser():
  description = "sgRSEA main stat function"
  epilog = "For command line options of each command, type %(prog)% COMMAND -h"
  argparser = ap.ArgumentParser(description=description, epilog = epilog)
  argparser.add_argument("-i","--input",dest = "infile",type=str,required=True, help = "sgRSEA input file, 4 columns")
  argparser.add_argument("-o","--output",dest = "outfile",type=str,required=True, help = "output file name")
  argparser.add_argument("--multiplier",dest = "multiplier",type=int, default = 50, help = "Multiplier to generate background")
  #argparser.add_argument("--bgtag",dest = "bgtag",type=str, default = "",help = "Sting to identify control sgRNAs")
  #argparser.add_argument("--bg-row-start",dest = "bgrowstart",type=int,default = -1, help = "Row count of the start of control sgRNA block")
  #argparser.add_argument("--bg-row-stop",dest = "bgrowstop",type=int, default=-1, help = "Row count of the stop of control sgRNA block")
  return(argparser)

#DEP#def getBackground(infile,nontag="",tagStart=0,tagStop=0):
#DEP#  '''Get background data frame from either string tag or row range.
#DEP#  '''
#DEP#  if len(nontag)>0:
#DEP#    bgFile = infile[infile.iloc[:,0].str.contains(nontag)]
#DEP#    dataFile = infile[infile.iloc[:,0].str.contains(nontag)==False]
#DEP#  elif tagStart<tagStop and tagStart<infile.shape[0]: #use the row range to get file
#DEP#    tagStart -= 1
#DEP#    tagStop = min(infile.shape[0],tagStop)
#DEP#    bgFile = infile.iloc[tagStart:tagStop,:]
#DEP#    dataFile = infile.iloc[0:tagStart-1,:].append(infile.iloc[tagStop:,:])
#DEP#  else:#Non background
#DEP#    logging.info("No designed NonTargeting sgRNAs, use input data as backgound.")
#DEP#    bgFile = copy.copy(infile)
#DEP#    dataFile = infile
#DEP#
#DEP#  logging.debug("bgFile: "+str(bgFile.shape[0]))
#DEP#  logging.debug("dataFile: "+str(dataFile.shape[0]))
#DEP#  if bgFile.shape[0]==0: #backgound file is empty
#DEP#    logging.info("There is no NonTargeting sgRNA found. Use input data as background.")
#DEP#    bgFile = copy.copy(infile)
#DEP#    dataFile = infile
#DEP#  return (dataFile,bgFile)

def addZstat(data_df, pNull):
  data_df['pmme'] = data_df['treat']/(data_df['treat']+data_df['ctrl'])
  data_df['zstat'] = data_df._get_numeric_data().apply(lambda x: zStat(x,pNull),axis=1)
  return data_df

def dataFilter(dataFile,sg_min = 1):
  '''Get table of sgRNA number per gene distribution. Filter out the genes with sgRNA number less than sg_min'''
  filter_data = dataFile.groupby('Gene').filter(lambda x : len(x)>=sg_min)
  #logging.debug(type(filter_data))
  #logging.debug(filter_data.columns.values)
  #filter_data['zstat'] = filter_data.loc[:,['treat_1','ctrl_1']].apply(lambda row: zStat(row.to_frame(['treat_1','ctrl_1']),0.4))
  gene_sample = filter_data.loc[:,'Gene']
  return (filter_data,gene_sample)


def pMME(treat, ctrl):
  '''1st numeric col is treatment and 2nd numeric col is control'''
  t_total = sum(treat)
  c_total = sum(ctrl)
  try:
    pmme = t_total * 1.0 / (t_total+c_total)
  except ZeroDivisionError:
    #logging.warning("ZeroDivisionError in pMME calculation, which indicates no reads count for treatment. Return 0 instead. Continuing")
    pmme = -1
  if not (0.0<=pmme<=1.0):
    logging.error("pMME range error. Exit")
    sys.exit(1)
  return pmme

def zStat(sg_row,pNull):
  #logging(sg_row)
  if pNull == 0 or pNull ==1:
    return None
  else:
    zscore = (sg_row.pmme-pNull) * math.sqrt(sg_row.treat+sg_row.ctrl) / math.sqrt(pNull*(1-pNull)) 
    return zscore

def maxMean(gene_zstat):
  scores = gene_zstat#.tolist()
  if len(scores)==0:
    logging.warning("No value for maxmean calculation")
    return 0
  sc_po = float(sum(scores[scores>0]))
  sc_ne = float(sum(scores[scores<0]))
  if abs(sc_po) >= abs(sc_ne):
    return sc_po/len(scores)
  else:
    return sc_ne/len(scores)

def getMatrixMaxmean(filtered_zdf):
  '''
     This function returns a dataframe with 'Gene', 'maxmean', 'sgCount'    
  '''
  maxmean_df = filtered_zdf.groupby("Gene").agg({'zstat':maxMean,'sgRNA':'count'}).reset_index()
  maxmean_df.columns = ["Gene","maxmean","sgcount"]
  #logging.debug(maxmean_df.head(10))
  return maxmean_df

#DEP#def sampleNull(genelist,bgFile):
#DEP#  sample_size = len(genelist)
#DEP#  samples = bgFile.iloc[np.random.randint(0,len(bgFile),size=sample_size)]
#DEP#  #make samples into gene_matrix structure list
#DEP#  samples['gene'] = genelist.values
#DEP#  logging.debug(samples.head(10))
#DEP#  #group_sample = samples.groupby('gene')
#DEP#  return samples

def splitPoints(census):
  '''Input a array with tuples (value, freq)'''
  split_at = [0]
  zscore_split_at = [0]
  value_list = []
  #census = [(2,4),(3,2),(6,1)]
  for value,freq in census:
    zscore_split_at.append(freq+zscore_split_at[-1])
    value_list.append(value)
    for i in range(freq):
      split_at.append(split_at[-1]+value)
  split_at.pop(0)
  split_at.pop(-1)
  zscore_split_at.pop(0)
  zscore_split_at.pop(-1)
  return (split_at, zscore_split_at, value_list)

def getMatrixMaxmean_null_numpy(datafile, multiplier=10, randSeed=None):
  #Get the partition of gene list
  zscorelist = np.array(datafile['zstat'])
  genecount = datafile['Gene'].value_counts()
  genebincount = np.bincount(genecount.values)
  genecount_census = zip(np.nonzero(genebincount)[0],genebincount[np.nonzero(genebincount)[0]])
  print genecount_census
  (split_at, zscore_split_at, sg_value_list) = splitPoints(genecount_census)
  print zscore_split_at
  np.random.seed(8512)
  logging.debug("start maxmean "+datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S'))
  maxmean_dic = {}
  for loop in range(multiplier):
    split_zscore = np.split(np.random.permutation(zscorelist),split_at)
    print len(split_zscore)
    zscore_group = np.split(split_zscore, zscore_split_at)
    #for m in range(1):
    #  print m,zscore_group[m]
    if len(zscore_group)!= len(sg_value_list):
      logging.error("Zscore group length does not match sgRNA group length. Exit")
      sys.exit(0)
    for i in range(len(sg_value_list)):
      k = sg_value_list[i]
      if not maxmean_dic.has_key(k):
        maxmean_dic[k]=[]
      #print type(zscore_group[i])
      #print zscore_group[i].tolist()
      #print np.asarray(zscore_group[i].tolist())
      
      result = (np.apply_along_axis(maxMean,1,zscore_group[i].tolist())).tolist()
      maxmean_dic[k]+=result
      #except:
      #  print "cannot maxmean",k
      #print maxmean_dic[2][0:10]
  logging.debug("stop maxmean "+datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S'))
#  for k,v in maxmean_dic.items():
#    print k,v[0:10]

def maxmeanSampleNull(genelist,bgFile,multiplier=10):
  '''Call sampleNull() multiplier times. Only keep maxmean dataframe to save memory.'''
  null_maxmean = pd.DataFrame()
  logging.debug("start copy df "+datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S'))
  null_df = copy.copy(bgFile)
  logging.debug("finish copy df "+ datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S'))
  np.random.seed(8512)
  nulls = []
  for i in range(multiplier):
    logging.debug("round "+str(i))
    logging.debug("start permutation:"+datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S'))
    newgene = np.random.permutation(genelist)
    null_df['Gene'] = newgene
    logging.debug("start calculate maxmean"+datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S'))
    newdf = getMatrixMaxmean(null_df)
    nulls.append(newdf)
  logging.debug("start concat df"+datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S'))
  null_maxmean = pd.concat(nulls)
  logging.debug("concat finished"+datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S'))
  #logging.debug("Null_maxmean_df")
  #logging.debug(null_maxmean[null_maxmean['sgcount'].isin([8,9])])
  return null_maxmean

def standardizeFactor(null_maxmean):
  '''Input: null_maxmean_df, Returns a df with 'sgcount', 'mean','std' '''
  group_null_maxmean = null_maxmean.groupby('sgcount')
  factor_df = pd.DataFrame({'mean':group_null_maxmean.maxmean.agg(np.mean),'std':group_null_maxmean.maxmean.agg(np.std)}).reset_index()
  return factor_df
  
#DEP#def standardize(datarow,factors):
#DEP#  '''Datarow is a pd.Series, index is ['gene','maxmean','sgcount']'''
#DEP#  factor = factors[factors['sgcount']==datarow.sgcount]
#DEP#  if factor.shape[0]==1:#Should be only one record
#DEP#    f_mean = factor.iloc[0]['mean']
#DEP#    f_std = factor.iloc[0]['std']
#DEP#  elif factor.shape[0]==0:
#DEP#    logging.warning("There is no record for genes with "+str(datarow.sgcount)+" sgRNAs. ")
#DEP#  else:#more than 1 record
#DEP#    logging.warning("There are more than 1 record for sgRNA count "+str(datarow.sgcount+" group. Check data."))
#DEP#  smm = (datarow.maxmean-f_mean)*1.0/f_std
#DEP#  return smm

def standardizeDF(maxmean_df,sFactors):
  #BC#maxmean_df['NScore'] = maxmean_df.apply(lambda x: standardize(x,sFactors),axis=1)
  df = maxmean_df.merge(sFactors,on="sgcount",how="left")
  df['NScore'] = (df['maxmean']-df['mean'])/df['std']
  return df
 #BC# return maxmean_df
  
def pvalue(tn,smm_null):
  '''smm_null is a dataframe'''
  pos_p = (smm_null[smm_null['NScore']>tn].shape[0]+1.0)/(1.0+smm_null.shape[0])
  return pos_p

def getPQ(data_maxmean_std,null_maxmean_std):
  data_maxmean_std['pos_p'] = data_maxmean_std.apply(lambda x: pvalue(x['NScore'],null_maxmean_std),axis=1)
  data_maxmean_std['neg_p'] = 1.0 -data_maxmean_std['pos_p'] 
  data_maxmean_std['pos_fdr'] = multipletests(data_maxmean_std.loc[:,'pos_p'],method='fdr_bh')[1]
  data_maxmean_std['neg_fdr'] = multipletests(data_maxmean_std.loc[:,'neg_p'],method='fdr_bh')[1]
  #data_maxmean_std['pos_rank'] = data_maxmean_std['pos_p'].rank().astype(int)
  #data_maxmean_std['neg_rank'] = data_maxmean_std['neg_p'].rank().astype(int)
  data_maxmean_std = data_maxmean_std.sort_values(by=['pos_fdr','NScore'],ascending=[True,False])
  data_maxmean_std = data_maxmean_std.reset_index(drop=True)
  data_maxmean_std['pos_rank'] = data_maxmean_std.index + 1
  data_maxmean_std['neg_rank'] = data_maxmean_std.shape[0] - data_maxmean_std['pos_rank']
  return data_maxmean_std

def runStatinfer(infile,outfile,multiplier):
  #if (len(nontag)!=0) or (tagStart>0 and tagStop>0):
  #  (dataFile,bgFile) = getBackground(infile,nontag,tagStart,tagStop)
  #  p0 = pMME(bgFile)
  #else:#Use dataset as background
  dataFile = pd.read_table(infile)
  p0 = pMME(dataFile.iloc[:,2],dataFile.iloc[:,3])
  if p0 ==0 or p0 ==1:
    logging.error("pMME for background equals to 0/1, indicating no counts for treatment or control. Please check your data. Exit.")
    sys.exit(1)
  #logging.debug(p0) 
  dataFile.columns = ['sgRNA','Gene','treat','ctrl']
  dataFile = addZstat(dataFile,p0)
  #logging.debug(dataFile.head(10))
  dataFile.to_csv("test_real_zscore.txt",sep="\t",header=True,index=False)
  (filtered_data,genelist) = dataFilter(dataFile,1)
  #dataStat = treat_group.size()
  logging.info("Calculating treatment maxmean...")
  data_maxmean_df = getMatrixMaxmean(filtered_data)
  data_maxmean_df.to_csv("test_real_maxmean.txt",sep="\t",header=True,index=False)
  logging.info("Sampling null distribution...")
#BC#  if isinstance(bgFile, pd.DataFrame): 
#BC#    bgFile = addZstat(bgFile, p0)
#BC#  else:
#BC#    bgFile = filtered_data
  null_maxmean_df = maxmeanSampleNull(genelist,filtered_data,multiplier)
  null_maxmean_df.to_csv("null_maxmean_df.xls",sep="\t",index=False)
  #logging.debug("Get standardize factors")
  factor_df = standardizeFactor(null_maxmean_df)
  #logging.debug(factor_df)
  logging.info("Standardization...")
  data_sdf = standardizeDF(data_maxmean_df,factor_df)
  null_sdf = standardizeDF(null_maxmean_df,factor_df)
  data_sdf.to_csv("test_real_standardize.txt",sep="\t",header=True,index=False)
  fdf = getPQ(data_sdf,null_sdf)
  fdf = fdf.loc[:,['Gene','sgcount','NScore','pos_p','pos_fdr','neg_p','neg_fdr','pos_rank','neg_rank']]
  fdf.to_csv(outfile,sep="\t",index=False)

def main():
  argparser = prepare_argparser()
  args = argparser.parse_args()
  infile = pd.read_table(args.infile)
  #print pMME(infile)
  df = getMatrixMaxmean_null_numpy(infile)
  #df.to_csv(args.outfile,sep="\t",index=False)
  #runStatinfer(args.infile,args.outfile,args.bgtag, args.bgrowstart, args.bgrowstop, args.multiplier)
  #runStatinfer(args.infile,args.outfile,args.multiplier)
  
if __name__ == '__main__':
	main()
	
