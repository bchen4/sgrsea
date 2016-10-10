#!/usr/bin/python
# programmer : bbc
# usage:

import sys
import os
import logging
import numpy as np
import pandas as pd
import argparse as ap
import sgcount
import normalization
import stattest
import reformatCountTable 
logging.basicConfig(level=10)

def prepare_argparser():
  description = "Run the whole sgRSEA suite"
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
  argparser.add_argument("--normalize_method",dest="method", default="total", type=str,help ="design file", choices=['total','median','upperquantile'])
  argparser.add_argument("--split-lib",dest = "splitlib",action='store_true', help="Lib A and B are sequenced separately")
  argparser.add_argument("-t","--treatment",dest="treat",type=str,required=True, help = "columns/name of treatment samples")
  argparser.add_argument("-c","--control",dest="ctrl",type=str,required=True, help="columns/name of control samples")
  argparser.add_argument("--t-lable",dest="treatlabel",type=str, help = "label of treatment samples")
  argparser.add_argument("--c-label",dest="ctrllabel",type=str, help="label of control samples")
  argparser.add_argument("--multiplier",dest = "multiplier",type=int, default = 30,required=True, help = "Multiplier to generate background")
  argparser.add_argument("--collapse-replicates",dest="collapsemethod",type=str, help = "Way to collapse replicates", default="auto", choices=['auto','stack','mean'])

  #argparser.add_argument("--bgtag",dest = "bgtag",type=str, default = "",help = "Sting to identify control sgRNAs")
  #argparser.add_argument("--bg-row-start",dest = "bgrowstart",type=int,default = -1, help = "Row count of the start of control sgRNA block")
  #argparser.add_argument("--bg-row-stop",dest = "bgrowstop",type=int, default=-1, help = "Row count of the stop of control sgRNA block")
  return argparser

def run(args):
  if args.infile != None:
    logging.info("This is one-file mode, only count and do noramlization")
  sgcount.run(agrs)
  if os.path.exists(args.outfile+".count.txt"):#count matrix file exsists
    normalization.normalization(args.outfile+".count.txt",args.outfile+".norm.txt",args.method, args.splitlib)
  if os.path.exists(args.outfile+".norm.txt"):
    files = reformatCountTable.runReformat(args.outfile+".norm.txt",args.designfile, args.outfile,args.treat, args.ctrl, args.collapsemethod)
  for fn in files:
    logging.info("Running test on "+fn)
    stattest.runStatinfer(fn,fn+".sgRSEA.xls",args.multiplier)

def main():
  argparser = prepare_argparser()
  args = argparser.parse_args()
  run(args)
  #normalization(args.infile, args.outfile, args.method, args.splitlib)
  #infile = pd.read_table(args.infile)
  #logging.debug("read file")
  #print norm(infile,"total")

if __name__ == '__main__':
  main()
  
