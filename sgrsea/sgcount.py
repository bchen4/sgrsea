# usage: count sgRNA from fastq file

import sys
from Bio import SeqIO
import argparse as ap
import logging
import pandas as pd
import numpy as np
from multiprocessing import Process, Queue
from multiprocessing import Pool
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
    lib_df = pd.DataFrame()
    (result, total_fq_count) = multicount_limit(dfile)
    (result_df, summary_df) = generatefinaltable(result, total_fq_count, lib, dfile)
  elif args.infile:
    #infile = open(args.infile,"r")
    lib = makelib([args.libfile],['sublib'])
    result = {}
    ###logging.debug("There are "+str(lib.shape[0])+" sgRNAs in the library")
    total_fq_count = {}
    (result_addr,total_count,sublib) = sgcount(args.infile, args.sgstart, args.sgstop, args.trim3,"count","sublib")
    result[sublib] = pd.read_table(result_addr+".tmpcount")
    total_fq_count[result_addr] = total_count
  #Get total reads count df
    (result_df, summary_df) = generatefinaltable(result, total_fq_count, lib, None)
  result_df.to_csv(args.outfile+".count.txt",sep="\t",index=False)
  summary_df.to_csv(args.outfile+".summary.txt",sep="\t",index=False)

def generatefinaltable(resultdic, totaldic, lib, dfile):
  '''
    merge df based on sublib  
  ''' 
  #Generate count table
  count_df = pd.DataFrame()
  for sublibname, c_df in resultdic.items():
    sublib = lib[sublibname]
    df = sublib.merge(c_df,on='Sequence',how='left')
    df = df.fillna(0)
    count_df = count_df.append(df)
  #Generate summary table
  count_columns = count_df.columns.tolist()
  count_columns.insert(0,count_columns.pop(count_columns.index('sgRNA')))
  count_columns.insert(1,count_columns.pop(count_columns.index('Gene')))
  count_columns.insert(2,count_columns.pop(count_columns.index('Sequence')))
  count_columns.insert(3,count_columns.pop(count_columns.index('sublib')))
  count_df = count_df.loc[:,count_columns]
  mapped_total = count_df.iloc[:,3:].groupby("sublib").sum().reset_index()
  mapped_total_df = pd.melt(mapped_total, id_vars=['sublib'],var_name=['sample'],value_name='mapped_reads')
  totalread_df = pd.DataFrame(totaldic.items(),columns=["filepath","total_reads"])
  if isinstance(dfile,pd.DataFrame):
    summary_df = dfile.merge(totalread_df,on="filepath")
    summary_df = summary_df.merge(mapped_total_df,on=['sublib','sample'])
  else:#single file
    summary_df = totalread_df
    summary_df = summary_df.join(mapped_total_df)
  summary_df['mapping_ratio'] = summary_df['mapped_reads']/summary_df['total_reads']
  summary_df = summary_df.loc[:,['filepath','sample','sublib','total_reads','mapped_reads','mapping_ratio']]
  return (count_df, summary_df)

def sgcountwrapper(args):
  #return fq file name, total count and sublit
  return sgcount(*args)

def multicount(dfile, number_of_workers=5):
  '''
  Parse design file, and call sgcount for each pair of files. 
  Combine all results to a Pandas dataframe and return,group by sublib
  '''
  work_num = number_of_workers
  work_pool = Pool(work_num)
  workers = []
  #make a list of arguments, each row/item is a call
  arglist = zip(dfile['filepath'].tolist(),dfile['sgstart'].tolist(),
      dfile['sgstop'].tolist(), dfile['trim3'].tolist(), 
      dfile['sample'].tolist(), dfile['sublib'].tolist())
  results = work_pool.map(sgcountwrapper, arglist)  
  #print results 
  #get results
  #logging.info("Retrive results")
  total_fq_count = {}
  samplefile = {}
  for info in results:
    (fqfilename,total_count,fsublib) = info#output should be a Pandas df
    total_fq_count[fqfilename] = total_count
    if not samplefile.has_key(fsublib):
      samplefile[fsublib] = []
    samplefile[fsublib].append(fqfilename)
  result_dic ={}
  for sublibname, files in samplefile.items():
    result = pd.read_table(files[0]+".tmpcount")
    try:
      os.remove(files[0]+".tmpcount")
    except:
      logging.warning("Cannot delete "+files[0]+".tmpcount")
    for fn in files[1:]:
      result = result.merge(pd.read_table(fn+".tmpcount"),on="Sequence",how="outer")
      try:
        os.remove(fn+".tmpcount")
      except:
        logging.warning("Cannot delete "+fqfilename+".tmpcount")
    result_dic[sublibname] = result
  return (result_dic, total_fq_count)

def makelib(libs, sublib):
  '''
  #Pandas dataframe columns are: ['sgRNA','Gene','Sequence','sublib']
  return a dictionary with sample as the key and df as value
  '''
  lib_dic = {}
  for i in range(len(libs)):
    df = pd.read_table(libs[i])
    df['sublib'] = [sublib[i]]*df.shape[0]
    seq_count = df['Sequence'].value_counts()
    seq_count_df = pd.DataFrame({'Sequence':seq_count.index, 'Count':seq_count.values})
    df_uniq = df[df['Sequence'].isin(seq_count_df[seq_count_df['Count']==1]['Sequence'])]
    lib_dic[sublib[i]] = df_uniq
  return lib_dic
   
def trimseq(seq,start,stop,trim3=None):
  if len(seq)> start > 0:
    new_start = start - 1
  else:
    new_start = 0
  if 0 < stop <= len(seq):
    new_stop  = stop
  else:
    new_stop = len(seq)
  if trim3:#find the pattern from tail and trim off the rest of seq
    trim_index = seq.rfind(trim3)
    if trim_index >0:
      new_stop = min(new_stop, trim_index+1)
  trim_seq = seq[new_start:new_stop]
  return trim_seq


def sgcount(fqfile,sgstart, sgstop, trim3, sample="count",sublib="sublib"):
  '''
  Count sequence frequency in a fastq file.
  Write result in a temp file due to Python multiprocessing hang with huge result. (It can only handle int and string, not a instance of a class)
  '''
  trim_out = open(fqfile+".trim","w")
  total_count = 0
  seqdic = {}
  for record in SeqIO.parse(fqfile,"fastq"):
    total_count += 1
    if total_count % 500000 ==0:
      logging.info("Processed "+fqfile+" "+str(total_count)+" reads...")
      break
    sequence =  trimseq(str(record.seq),sgstart, sgstop, trim3)
    print >> trim_out, ">"+record.id+"\n"+sequence
    if not seqdic.has_key(sequence):
      seqdic[sequence] = 0
    seqdic[sequence]+=1
  seq_count = pd.DataFrame(seqdic.items(),columns=['Sequence',sample])
  #logging.info(seq_count.shape)
  seq_count.to_csv(fqfile+".tmpcount",sep="\t",index=False)
  trim_out.close()
  return (fqfile,total_count,sublib)

def run(args):
  runsgcount(args)

def main():
  argparser = prepare_argparser()
  args = argparser.parse_args()
  runsgcount(args)

if __name__=="__main__":
  main()
