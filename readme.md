#sgRSEA: Enrichment Analysis of CRISPR/Cas9 Knockout Screen Data

sgRSEA is a tool to identify significant genes in GeCKO experiments.

Version:alpha

Last modified: 10/10/2016

Authors: Beibei Chen (beibei.chen@utsouthwestern.edu), Jungsik Noh (junsik.noh@utsouthwestern.edu)

Maintainer: Beibei Chen

[Download the latest stable version of sgRSEA] (http://)
[Download the R package of sgRSEA] (https://cran.r-project.org/web/packages/sgRSEA/index.html)

For detailed [manual](https://github.com/bchen4/sgrsea/wiki/Manual)

##Prerequisites
* Python 2.7
* [Pandas](http://pandas.pydata.org/)
* [Numpy] (http://www.numpy.org/)
* stattest 

##Installation 

###Option 1: from source code
1 Download sgRSEA
```bash
git clone
cd sgRSEA
```

2 Install required packages
```bash
pip install -r requirements.txt
```

3 Install sgRSEA
```python
python setup.py install
```

##Usage

sgRSEA can be run on specific steps as well as fastq-to-result using a single command.

```
sgRSEA: identify significant genes in CRISPR-Cas9 experiment

positional arguments:
  {count,reformat,run,stattest,normalization}
    run                 Run the whole program from fastq to result
    count               Get sgRNA count matrix
    normalization       Normalize sgRNA count matrix
    reformat            Reformat count table for sgRSEA stat test
    stattest            Identify significant genes from normalized count table

optional arguments:
  -h, --help            show this help message and exit

For command line options of each command, type sgrsea COMMAND -h
```

###Get count matrix from fastq files
**sgRSEA count** can be used on single fastq file and multiple fastq files. When using multiple fastq files, a [design file](http://) containing design information should be provided.

```bash
usage: sgrsea count [-h] [-i INFILE | -d DESIGNFILE] -o OUTFILE [-l LIBFILE]
                    [--sgstart SGSTART] [--sgstop SGSTOP] [--trim3 TRIM3]

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --input INFILE
                        input fastq
  -d DESIGNFILE, --design DESIGNFILE
                        design file
  -o OUTFILE, --output OUTFILE
                        output
  -l LIBFILE, --library LIBFILE
                        Gene locus in bed format
  --sgstart SGSTART     The first nucleotide sgRNA starts. 1-index
  --sgstop SGSTOP       The last nucleotide sgRNA starts. 1-index
  --trim3 TRIM3         The trimming pattern from 3'. This pattern and the
                        following sequence will be removed
```

###Normalize sgRNA count matrix

```bash
usage: sgrsea normalization [-h] -i INFILE
                            [--normalize-method {total,median,upperquantile}]
                            -o OUTFILE [--split-lib]

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --input INFILE
                        input count file matrix
  --normalize-method {total,median,upperquantile}
                        design file
  -o OUTFILE, --output OUTFILE
                        output
  --split-lib           Lib A and B are sequenced separately
```

###Convert data matrix to sgRSEA input

sgRSEA takes a data matrix with 4 columns: sgRNA, Gene, treatment, control.
**sgreas reformat** to collapse replicates and make multiple input files for *stattest* if there are more than 1 comparison needed to be done accroding to the design file.

```bash
usage: sgrsea reformat [-h] -i INFILE -o OUTFILE [-d DESIGNFILE] -t TREAT -c
                       CTRL [--collapse-replicates {auto,stack,mean}]

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --input INFILE
                        input BAM file
  -o OUTFILE, --output OUTFILE
                        output
  -d DESIGNFILE, --design DESIGNFILE
                        output
  -t TREAT, --treatment TREAT
                        columns/name of treatment samples
  -c CTRL, --control CTRL
                        columns/name of control samples
  --collapse-replicates {auto,stack,mean}
                        Way to collapse replicates
```

###Find significant genes for each comparison

```bash
usage: sgrsea stattest [-h] -i INFILE -o OUTFILE --multiplier MULTIPLIER

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --input INFILE
                        sgRSEA input file, 4 columns
  -o OUTFILE, --output OUTFILE
                        output file name
  --multiplier MULTIPLIER
                        Multiplier to generate background
```

##Examples

Run the suite:
```python
sgrsea run -d ../test/fulldata/design.txt -o ../test/fulldata/testfull --multiplier 100 -t BR,DI -c UN,UN --split-lib
```


##Input

###[Fastq file](https://en.wikipedia.org/wiki/FASTQ_format)

###Library file
Library file should have 3 columns in the following order, with the column names exactly the same: 
* Gene 
* sgRNA
* Sequence

Example:
```
Gene  sgRNA             Sequence
A1BG  HGLibA_00001  GTCGCTGAGCTCCGATTCGA
A1BG  HGLibA_00002  ACCTGTAGTTGCCGGCGTGC
A1BG  HGLibA_00003  CGTCAGCGTCACATTGGCCA
A1CF  HGLibA_00004  CGCGCACTGGTCCAGCGCAC
A1CF  HGLibA_00005  CCAAGCTATATCCTGTGCGC
A1CF  HGLibA_00006  AAGTTGCTTGATTGCATTCT
A2M   HGLibA_00007  CGCTTCTTAAATTCTTGGGT
A2M   HGLibA_00008  TCACAGCGAAGGCGACACAG
A2M   HGLibA_00009  CAAACTCCTTCATCCAAGTC
A2ML1 HGLibA_00010  AAATTTCCCCTCCGTTCAGA
```




##License
See the LICENSE file for license rights and limitations (MIT).















https://cran.r-project.org/web/packages/sgRSEA/index.html
