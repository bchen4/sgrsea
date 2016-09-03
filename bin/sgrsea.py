import sys
sys.path.insert(1,'../')

import argparse
import logging
import os
import pandas as pd
import numpy as np
import locale

from crispr import *
logging.basicConfig(format='%(levelname)s:%(message)s',level=logging.DEBUG)


def prepare_argparser():
  description = "Look for enriched genes in CRISPR-Cas9 experiment"
  epilog = "For command line options of each command, type %(prog)s COMMAND -h"
  argparser = argparse.ArgumentParser(description=description, epilog = epilog)
  argparser.add_argument("-i","--input",dest = "infile", type = str, required = True, help = "Input count table")
  argparser.add_argument("-t","--treatment",dest = "treatcol", type = str, required = True, default = 3 ,help = "The column numbers of treatment. -t 3,4,5")
  argparser.add_argument("-c","--control",dest="ctrlcol",type=str,required=True, default=4, help = "The column number of control. -c 6,7,8")
  argparser.add_argument("-o","--output",dest = "outfile", type = str,required = False, help = "output file, default is a input.crispr.xls file in current folder")
  argparser.add_argument("-n","--normalization",dest = "normalizeMethod", type = str,required = True, default="total", help = "normalization method: upperquartile or total, default is normalize the data by total reads count")
  argparser.add_argument("-b","--nontarget",dest = "nonT", type = str ,required = False, help = "The label to identify nontarget sgRNAs. Exlusive to -r. You must set either -b or -r.")
  argparser.add_argument("-r","--nontargetRow",dest = "nonTrow", type = str,required = False, help = "Row ranges of nontargeting sgRNAs. Format: start,end. If only one number is provided, that number will be used as starting row and file end will be assumed as the last row for nontargeting sgRNAs")
  argparser.add_argument("-m","--multiplier",dest="multiplier", type = int, default = 10, help = "Count of input file sized Null distribution. Default = 10")
  argparser.add_argument("-p","--pvalue",dest = "pvalue", type = float,required = False, default = 0.05, help = "FDR cutoff for significant genes. Default is 0.05")
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
	#Validating parameters
  nontag = ""
  start = 0
  stop = 0
  if args.nonT != None:
    nontag = args.nonT
  else:
    if args.nonTRow == None:
      logging.error("You must provide either -b or -r. Exiting program.")
      sys.exit(0)
    else:#nontargting row range provided
      rows = args.nonTRow.split(",")
      if len(rows)==2:
        start = int(rows[0])
        stop = int(rows[1])
      else:
        if start > stop:
          logging.error("NonTarget row start is bigger than stop, please check your parameters. Exiting. ")
          sys.exit(0)
        else:
          start = max(0,int(rows[0]))
          stop = 10**20 #set to a big number so the end of file will be used
  #Validating treatments and controls tags 	
  treatments = args.treatcol.split(",")
  controls = args.ctrlcol.split(",")
  if len(treatments) ==0 or len(controls)==0:
    logging.error("Please provide treatments and controls column numbers. Exit.")
    sys.exit(1)
  elif len(np.intersect1d(treatments,controls))>0:
    logging.error("Treatments and controls column numbers overlap. Please check your file. Exit.")
    sys.exit(1)
  elif min(treatments) < 3 or min(controls) < 3:
    logging.error("Count data should start from the 3rd column. Columns received contain numbers smaller than 3. Please check your file. Exit.")
    sys.exit(1)
  
  logging.debug("Finish check nonT info")
  runCrispr(args,treatments,controls,nontag,start,stop)

if __name__=="__main__":
	main()
