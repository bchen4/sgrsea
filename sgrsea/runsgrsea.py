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
from multiprocessing import Process, Queue
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
  argparser.add_argument("--num-threads",dest="threads",type=int,default = 1,help = "Number of threads to use.")
  argparser.add_argument("--normalize_method",dest="method", default="total", type=str,help ="design file", choices=['total','median','upperquantile'])
  argparser.add_argument("--split-lib",dest = "splitlib",action='store_true', help="Lib A and B are sequenced separately")
  argparser.add_argument("-t","--treatment",dest="treat",type=str,required=True, help = "columns/name of treatment samples")
  argparser.add_argument("-c","--control",dest="ctrl",type=str,required=True, help="columns/name of control samples")
  argparser.add_argument("--t-lable",dest="treatlabel",type=str, help = "label of treatment samples")
  argparser.add_argument("--c-label",dest="ctrllabel",type=str, help="label of control samples")
  argparser.add_argument("--multiplier",dest = "multiplier",type=int, default = 50, help = "Multiplier to generate background")
  argparser.add_argument("--random-seed",dest = "randomSeed",type=int, default = None, help = "Random seed to control permutation process")
  argparser.add_argument("--collapse-replicates",dest="collapsemethod",type=str, help = "Way to collapse replicates", default="auto", choices=['auto','stack','mean'])
  argparser.add_argument("--no-count", dest="nocount", default=False, action='store_true', help="Skip counting step. Uses output contents as input")

  #argparser.add_argument("--bgtag",dest = "bgtag",type=str, default = "",help = "Sting to identify control sgRNAs")
  #argparser.add_argument("--bg-row-start",dest = "bgrowstart",type=int,default = -1, help = "Row count of the start of control sgRNA block")
  #argparser.add_argument("--bg-row-stop",dest = "bgrowstop",type=int, default=-1, help = "Row count of the stop of control sgRNA block")
  return argparser

def run(args):
  if isinstance(args.infile, str):
    logging.info("This is one-file mode, only count and do noramlization")
  #logging.info("Start to count fastq")
  if not args.nocount:
      sgcount.run(args)
  if os.path.exists(args.outfile+".count.txt"):#count matrix file exsists
    normalization.normalization(args.outfile+".count.txt",args.outfile+".norm.txt",args.method, args.splitlib)
  else:
    logging.error("There is no count file. Exit.")
    sys.exit(1)
  if isinstance(args.designfile,str):#multiple files, go further
    if os.path.exists(args.outfile+".norm.txt"):
      files = reformatCountTable.runReformat(args.outfile+".norm.txt",args.designfile, args.outfile,args.treat, args.ctrl, args.collapsemethod)
    else:
      logging.error("There is no normalized count file. Exit.")
      sys.exit(1)
    if len(files)==0:
      logging.error("There is no input files for sgRSEA to run. Exit.")
      sys.exit(1)
    else:
      work_num = len(files)
      work_queue = Queue()
      workers = []
      for fn in files:
        logging.info("Running test on "+fn)
        p = Process(target = callstat, args=(work_queue,fn,fn+".sgRSEA.xls",args.multiplier, args.randomSeed))
        workers.append(p)
        p.start()
      for process in workers:
        process.join()

def callstat(queue, fname, outname, multiplier, randomseed):
  queue.put(stattest.runStatinfer(fname,outname,multiplier,randomseed))
  
def main():
  argparser = prepare_argparser()
  args = argparser.parse_args()
  run(args)

if __name__ == '__main__':
  main()
  
