#!/usr/bin/python
# programmer : bbc
# usage:

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
    lib.to_csv(args.outfile,sep="\t",index=False)
    result = multicount(dfile)
  elif args.infile:
    infile = open(args.infile,"r")
    lib = makelib([args.libfile],['sublib'])
    #logging.debug("There are "+str(lib.shape[0])+" sgRNAs in the library")
    result = sgcount(infile, args.sgstart, args.sgstop, args.trim3)
  #merge results df with lib
  result_df = lib.merge(result,on="Sequence",how="left")
  result_df = result_df.fillna(0)
  result_df.to_csv(args.outfile+".count.txt",sep="\t",index=False)


def callsgcount(queue,fname,sgstart,sgstop,trim3,label):
  print "file name", fname
  queue.put(sgcount(fname,sgstart,sgstop,trim3,label))

def testfunc(f,a,b):
  return f+":"+str(a+b)

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

  for process in workers:
    process.join()

  result = work_queue.get()#output should be a Pandas df
  for i in range(work_num-1):
    result = result.merge(work_queue.get(),on="Sequence",how="outer")
  return result


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
  #logging.debug(seq)
  #logging.debug(start)
  #logging.debug(stop)
  #logging.debug(len(seq))
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
  mapped = 0
  total_count = 0
  seqdic = {}
  for record in SeqIO.parse(fqfile,"fastq"):
    total_count += 1
    sequence =  trimseq(str(record.seq),sgstart, sgstop, trim3)
    if not seqdic.has_key(sequence):
      seqdic[sequence] = 0
    seqdic[sequence]+=1
  seq_count = pd.DataFrame(seqdic.items(),columns=['Sequence',label])
  return seq_count

def main():
  argparser = prepare_argparser()
  args = argparser.parse_args()
  runsgcount(args)
  #makelib_pd(["demolib.txt"],["A"])
  #df = pd.read_table("demolib.txt")
  #makelib_pd(df)

if __name__=="__main__":
  main()
