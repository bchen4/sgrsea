import sys
sys.path.insert(1,'../')

import argparse
import logging
import os
import pandas as pd
import numpy as np
import locale
from sgrsea import *
from sgrsea import sgcount_fastq
import inspect
print inspect.getfile(sgrsea)

logging.basicConfig(format='%(levelname)s:%(message)s',level=logging.DEBUG)


def prepare_argparser():
  description = "sgRSEA: identify significant genes in CRISPR-Cas9 experiment"
  epilog = "For command line options of each command, type %(prog)s COMMAND -h"
  argparser = argparse.ArgumentParser(description=description, epilog = epilog)
#  argparser.add_argument("--version", action="version", version="%(prog)s"+sgrsea_VERSION)
  subparsers = argparser.add_subparsers( dest = 'subcommand_name' )
  subc_run = subparsers.add_parser("run", help="Run the whole program from fastq to result")
  subc_count = subparsers.add_parser("count", help="Get sgRNA count matrix")
  subc_normalize = subparsers.add_parser("normalize", help="Normalize count matrix")
  subc_stattest = subparsers.add_parser("stattest", help="Identify significant genes from normalized count table")
#BC#  argparser.add_argument("-i","--input",dest = "infile", type = str, required = True, help = "Input count table")
#BC#  argparser.add_argument("-t","--treatment",dest = "treatcol", type = str, required = True, default = 3 ,help = "The column numbers of treatment. -t 3,4,5")
#BC#  argparser.add_argument("-c","--control",dest="ctrlcol",type=str,required=True, default=4, help = "The column number of control. -c 6,7,8")
#BC#  argparser.add_argument("-o","--output",dest = "outfile", type = str,required = False, help = "output file, default is a input.crispr.xls file in current folder")
#BC#  argparser.add_argument("-n","--normalization",dest = "normalizeMethod", type = str,required = True, default="total", help = "normalization method: upperquartile or total, default is normalize the data by total reads count")
#BC#  argparser.add_argument("-b","--nontarget",dest = "nonT", type = str ,required = False, help = "The label to identify nontarget sgRNAs. Exlusive to -r. You must set either -b or -r.")
#BC#  argparser.add_argument("-r","--nontargetRow",dest = "nonTrow", type = str,required = False, help = "Row ranges of nontargeting sgRNAs. Format: start,end. If only one number is provided, that number will be used as starting row and file end will be assumed as the last row for nontargeting sgRNAs")
#BC#  argparser.add_argument("-m","--multiplier",dest="multiplier", type = int, default = 10, help = "Count of input file sized Null distribution. Default = 10")
#BC#  argparser.add_argument("-p","--pvalue",dest = "pvalue", type = float,required = False, default = 0.05, help = "FDR cutoff for significant genes. Default is 0.05")
  return(argparser)



def runCrispr(args,treatments,controls,tag="",tagStart=0,tagStop=0):
  try:
    infile = pd.read_table(args.infile,sep="\t")
  except:
    logging.error("Can not open input file, please check your file name. Exiting program.")
    sys.exit(1)
  if infile.shape[0]<=0:# No rows in input file
    logging.error("There is no data in the file. Exit.")
    sys.exit(0)
	#Normalization
  infile.fillna(0) #replace missing data with 0
  norm_flag = "single"
  if not (len(treatments)==1 and len(controls) ==1): # replicates exsit 
    norm_flag = "multiple"
  logging.debug(infile.columns.values) 
  normfile = Normalization.norm(infile,args.normalizeMethod,norm_flag)
  logging.debug("Normalization finished.")
  logging.debug(normfile.head(10))
	#Statistical test
  result = StatInfer.runTest(normfile,treatments,controls,tag,tagStart,tagStop, args.multiplier)
	#Output
  result.to_csv(args.outfile,sep="\t",header=True,index=False)



def main():
  logging.debug("main")
  arg_parser = prepare_argparser()
  args = arg_parser.parse_args()
  subcommand  = args.subcommand_name
  if subcommand == "count":
    import sgcount_fastq
    sgrsea.sgcount_fastq.run(args)
  elif subcommand == "normalize":
    sgrsea.Normalization.run(args)
  elif subcommand == "stattest":
    sgrsea.stattest.run(args)
  elif subcommand == "run":
    sgrsea.runsgrsea.run(args)
	#Validating parameters
#DEP#  nontag = ""
#DEP#  start = 0
#DEP#  stop = 0
#DEP#  if args.nonT != None:
#DEP#    nontag = args.nonT
#DEP#  else:
#DEP#    if args.nonTRow == None:
#DEP#      logging.error("You must provide either -b or -r. Exiting program.")
#DEP#      sys.exit(0)
#DEP#    else:#nontargting row range provided
#DEP#      rows = args.nonTRow.split(",")
#DEP#      if len(rows)==2:
#DEP#        start = int(rows[0])
#DEP#        stop = int(rows[1])
#DEP#      else:
#DEP#        if start > stop:
#DEP#          logging.error("NonTarget row start is bigger than stop, please check your parameters. Exiting. ")
#DEP#          sys.exit(0)
#DEP#        else:
#DEP#          start = max(0,int(rows[0]))
#DEP#          stop = 10**20 #set to a big number so the end of file will be used
#DEP#  #Validating treatments and controls tags 	
#DEP#  treatments = args.treatcol.split(",")
#DEP#  controls = args.ctrlcol.split(",")
#DEP#  if len(treatments) ==0 or len(controls)==0:
#DEP#    logging.error("Please provide treatments and controls column numbers. Exit.")
#DEP#    sys.exit(1)
#DEP#  #elif len(np.intersect1d(treatments,controls))>0:
#DEP#  #  logging.error("Treatments and controls column numbers overlap. Please check your file. Exit.")
#DEP#  #  sys.exit(1)
#DEP#  #elif min(treatments) < 3 or min(controls) < 3:
#DEP#  #  logging.error("Count data should start from the 3rd column. Columns received contain numbers smaller than 3. Please check your file. Exit.")
#DEP#    sys.exit(1)
#DEP#  
#DEP#  logging.debug("Finish check nonT info")
#DEP#  runCrispr(args,treatments,controls,nontag,start,stop)

if __name__=="__main__":
	main()
