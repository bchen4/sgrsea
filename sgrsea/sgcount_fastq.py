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
  #argparser.add_argument("-o","--output",dest = "outfile",type=str,required=True, help="output")
  #argparser.add_argument("-l","--library",dest="libfile",type=str,required=True, help = "Gene locus in bed format")
  #argparser.add_argument("--sgstart",dest="sgstart",type=int, default=-1,help = "The first nucleotide sgRNA starts. 1-index")
  #argparser.add_argument("--sgstop",dest="sgstop",type=int, default=-1,help = "The last nucleotide sgRNA starts. 1-index")
  #argparser.add_argument("--trim3",dest="trim3",type=str,help = "The trimming pattern from 3'. This pattern and the following sequence will be removed")
  return(argparser)

def runsgcount(args):
  if args.infile and args.designfile:
    logging.error("You must provide infile OR design file. Quit.")
    sys.exit(1)
  if args.designfile:#parse design file
    #If design file is provided, sgRNA start, stop, trim3 will be ignored
    dfile = pd.read_table(args.designfile)
    result = multicount(args)
  elif args.infile:
    infile = open(args.infile,"r")
    lib = makelib(args.libfile)
    #outfile = open(args.outfile+".count.txt","w")
    #outsummary = open(args.outfile+".count.summary","w")
    result = sgcount(infile, lib, args.sgstart, args.sgstop, args.trim3)

def testfunc(queue,pname,fname):
  print "Process name",pname,"file name", fname
  queue.put(pname)

def multicount(args):
  infile = open(args.designfile,"r")
  filebuf = infile.readlines()
  work_num = len(filebuf)
  work_queue = Queue()
  workers = []
  for index in range(work_num):
    p = Process(target=testfunc, args=(work_queue,index, filebuf[index]))
    workers.append(p)
    p.start()

  for process in workers:
    process.join()

  result = []
  for i in range(work_num):
    result.append(work_queue.get())
  print result

#DEL#def makelib(*libs):
#DEL#  libdic = {}
#DEL#  dupseq = []
#DEL#  for libfile in libs:
#DEL#    for row in libfile:
#DEL#      buf = row.rstrip().split("\t")
#DEL#      if not libdic.has_key(buf[1]):
#DEL#        libdic[buf[1]]=[buf[0],buf[2],0]
#DEL#      else:
#DEL#        dupseq.append(buf[1])
#DEL#    for k in dupseq:
#DEL#      libdic.pop(k,None)
#DEL#  return libdic

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
  return df_uniq.head()
    
def sgcount(fqfile, lib, sgstart, sgstop, trim3):
  mapped = 0
  total_count = 0
  for record in fqfile:
    total_count +=1
    #if i % 1000000 ==0:
  #    logging.info("Processed "+str(i% 1000000)+"M reads...")
    if libdic.has_key(str(record.seq)):
      mapped+=1
      libdic[str(record.seq)][2]+=1
  for values in libdic.values():
    print >> outfile,"%s\t%s\t%d\t%f" % (values[0],values[1],values[2],(values[2]*1000000.0/mapped)+1)
  print >>outsummary,"%s\t%d\t%d\t%f" % (args.outfile,total_count,mapped,mapped*1.0/total_count)
  outfile.close()


def main():
  argparser = prepare_argparser()
  args = argparser.parse_args()
  #runsgcount(args)
  makelib_pd(["demolib.txt"],["A"])
  #df = pd.read_table("demolib.txt")
  #makelib_pd(df)

if __name__=="__main__":
  main()
