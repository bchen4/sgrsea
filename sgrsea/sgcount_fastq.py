#!/usr/bin/python
# programmer : beibei chen
# usage: count sgRNA from fastq file
# last modification: 9/9/2016

import sys
from Bio import SeqIO
import argparse as ap
import logging
import pandas as pd
import numpy as np
from multiprocessing import Process, Queue
logging.basicConfig(level=10)


def prepare_argparser():
  description = "Get sgRNA counts from fastq file. For multiple files, use design file -d."
  epilog = "For command line options of each command, type %(prog)% COMMAND -h"
  argparser = ap.ArgumentParser(description=description, epilog = epilog)
  group = argparser.add_mutually_exclusive_group()
  group.add_argument("-i","--input",dest = "infile",type=str, help="input fastq")
  group.add_argument("-d","--design",dest="designfile",type=str,help = "design file")
  argparser.add_argument("-o","--output",dest = "outfile",type=str,required=True, help="output")
  argparser.add_argument("-l","--library",dest="libfile",type=str, help = "Gene locus in bed format")
  argparser.add_argument("--sgstart",dest="sgstart",type=int, default=-1,help = "The first nucleotide sgRNA starts. 1-index")
  argparser.add_argument("--sgstop",dest="sgstop",type=int, default=-1,help = "The last nucleotide sgRNA starts. 1-index")
  argparser.add_argument("--trim3",dest="trim3",type=str,help = "The trimming pattern from 3'. This pattern and the following sequence will be removed")
  return(argparser)

def runsgcount(args):
  if args.infile and args.designfile:
    logging.error("You must provide infile OR design file. Quit.")
    sys.exit(1)
  if args.designfile:#parse design file
    #If design file is provided, sgRNA start, stop, trim3 will be ignored
    dfile = pd.read_table(args.designfile)
    libinfo = dfile.loc[:,['lib','sublib']].drop_duplicates()
    lib = makelib(libinfo['lib'].tolist(),libinfo['sublib'].tolist())
    #lib.to_csv(args.outfile,sep="\t",index=False)
    (result, total_fq_count) = multicount(dfile)
  elif args.infile:
    infile = open(args.infile,"r")
    lib = makelib([args.libfile],['sublib'])
    ##logging.debug("There are "+str(lib.shape[0])+" sgRNAs in the library")
    total_fq_count = {}
    (result_addr,total_count) = sgcount(infile, args.sgstart, args.sgstop, args.trim3)
    result = pd.read_table(fqfile+".tmpcount")
    total_fq_count[fqfile] = total_count
  #Get total reads count df
  result_total_df = pd.DataFrame(total_fq_count.items(),columns=["filepath","total_reads"])
  summary_df = dfile.merge(result_total_df,on="filepath")
  #merge results df with lib
  result_df = lib.merge(result,on="Sequence",how="left")
  result_df = result_df.fillna(0)
  #logging.debug(result_df.head())
  mapped_total = result_df.iloc[:,3:].groupby("sublib").sum().reset_index()
  mapped_total_df = pd.melt(mapped_total, id_vars=['sublib'],var_name=['label'],value_name='mapped_reads')
  summary_df = summary_df.merge(mapped_total_df,on=['sublib','label'])
  summary_df['mapping_ratio'] = summary_df['mapped_reads']/summary_df['total_reads']
  summary_df = summary_df.loc[:,['filepath','label','sublib','total_reads','mapped_reads','mapping_ratio']]
  result_df.to_csv(args.outfile+".count.txt",sep="\t",index=False)
  summary_df.to_csv(args.outfile+".summary.txt",sep="\t",index=False)

def callsgcount(queue,fname,sgstart,sgstop,trim3,label):
  queue.put(sgcount(fname,sgstart,sgstop,trim3,label))

def multicount(dfile):
  '''
  Parse design file, and call sgcount for each file. 
  Combine all results to a Pandas dataframe and return
  '''
  work_num = dfile.shape[0]
  work_queue = Queue()
  workers = []
  for index in range(work_num):
    record = dfile.iloc[index,:]
    p = Process(target=callsgcount, args=(work_queue,record['filepath'],
      record['sgstart'],record['sgstop'],record['trim3'],record['label']))
    workers.append(p)
    p.start()
  
  logging.info("Start join process")
  for process in workers:
    process.join()
  logging.info("All process finished") 
  #get results
  logging.info("Retrive results")
  total_fq_count = {}
  (fqfilename,total_count) = work_queue.get()#output should be a Pandas df
  result = pd.read_table(fqfilename+".tmpcount")
  try:
    os.remove(fqfilename+".tmpcount")
  except:
    logging.warning("Cannot delete "+fqfilename+".tmpcount")
  total_fq_count[fqfilename] = total_count
  for i in range(work_num-1):
    (fqfilename,total_count) = work_queue.get()#output should be a Pandas df
    result = result.merge(pd.read_table(fqfilename+".tmpcount"),on="Sequence",how="outer")
    try:
      os.remove(fqfilename+".tmpcount")
    except:
      logging.warning("Cannot delete "+fqfilename+".tmpcount")
    total_fq_count[fqfilename] = total_count
  return (result, total_fq_count)


def makelib(libs, sublib):
  #Pandas dataframe columns are: ['sgRNA','Gene','Sequence','sublib']
  df = pd.DataFrame(columns=['sgRNA','Gene','Sequence','sublib'])
  for i in range(len(libs)):
    newdf = pd.read_table(libs[i])
    newdf['sublib'] = [sublib[i]]*newdf.shape[0]
    df = df.append(newdf)
  seq_count = df['Sequence'].value_counts()
  seq_count_df = pd.DataFrame({'Sequence':seq_count.index, 'Count':seq_count.values})
  df_uniq = df[df['Sequence'].isin(seq_count_df[seq_count_df['Count']==1]['Sequence'])]
  return df_uniq
   
def trimseq(seq,start,stop,trim3=None):
  if len(seq)> start > 0:
    new_start = start - 1
  if 0 < stop <= len(seq):
    new_stop  = stop
  if trim3:#find the pattern from tail and trim off the rest of seq
    trim_index = seq.rfind(trim3)
    if trim_index >0:
      new_stop = min(new_stop, trim_index+1)
  trim_seq = seq[new_start:new_stop]
  return trim_seq


def sgcount(fqfile,sgstart, sgstop, trim3,label="count"):
  '''
  Count sequence frequency in a fastq file.
  Write result in a temp file due to Python multiprocessing hang with huge result. (It can only handle int and string, not a instance of a class)
  '''
  total_count = 0
  seqdic = {}
  for record in SeqIO.parse(fqfile,"fastq"):
    total_count += 1
    if total_count % 500000 ==0:
      logging.info("Processed "+fqfile+" "+str(total_count)+" reads...")
      #break
    sequence =  trimseq(str(record.seq),sgstart, sgstop, trim3)
    if not seqdic.has_key(sequence):
      seqdic[sequence] = 0
    seqdic[sequence]+=1
  seq_count = pd.DataFrame(seqdic.items(),columns=['Sequence',label])
  logging.info(seq_count.shape)
  seq_count.to_csv(fqfile+".tmpcount",sep="\t",index=False)
  return (fqfile,total_count)

def main():
  argparser = prepare_argparser()
  args = argparser.parse_args()
  runsgcount(args)

if __name__=="__main__":
  main()
