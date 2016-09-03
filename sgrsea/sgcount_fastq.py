#!/usr/bin/python
# programmer : bbc
# usage:

import sys
from Bio import SeqIO
import argparse as ap
import logging
import pandas as pd
import numpy as np

logging.basicConfig(level=10)


def prepare_argparser():
  description = "Get sgRNA counts from fastq file"
  epilog = "For command line options of each command, type %(prog)% COMMAND -h"
  argparser = ap.ArgumentParser(description=description, epilog = epilog)
  argparser.add_argument("-i","--input",dest = "infile",type=str,required=True, help="design file")
  argparser.add_argument("-o","--output",dest = "outfile",type=str,required=True, help="output")
  argparser.add_argument("-l","--library",dest="libfile",type=str,required=True, help = "Gene locus in bed format")
  #argparser.add_argument("--sgstart",dest="sgstart",type=int, default=-1,help = "The first nucleotide sgRNA starts. 1-index")
  #argparser.add_argument("--sgstop",dest="sgstop",type=int, default=-1,help = "The last nucleotide sgRNA starts. 1-index")
  argparser.add_argument("--trim3",dest="trim3",type=str,help = "The trimming pattern from 3'. This pattern and the following sequence will be removed")
  return(argparser)


def main():
  argparser = prepare_argparser()
  args = argparser.parse_args()

  try:
    infile = SeqIO.parse(args.infile,"fastq")
    lib = open(args.libfile,"r")
    outfile = open(args.outfile+".count.txt","w")
    outsummary = open(args.outfile+".count.summary","w")
  except IOError,message:
    print >> sys.stderr, "cannot open file",message
    sys.exit(1)
  libdic = {}
  dupseq = []
  for row in lib:
    buf = row.rstrip().split("\t")
    if not libdic.has_key(buf[1]):
      libdic[buf[1]]=[buf[0],buf[2],0]
    else:
      dupseq.append(buf[1])
  for d in dupseq:
    print d
  for k in dupseq:
    libdic.pop(k,None)
  mapped = 0
  total_count = 0
  for record in infile:
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

if __name__=="__main__":
  main()
