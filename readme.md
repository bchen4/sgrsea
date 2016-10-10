#sgRSEA: identify significant genes in GeCKO

sgRSEA is a tool to identify significant genes in GeCKO experiments.

Version:alpha

Last modified: 10/10/2016

Authors: Beibei Chen (beibei.chen@utsouthwestern.edu), Jungsik Noh (junsik.noh@utsouthwestern.edu)

Maintainer: Beibei Chen

[Download the latest stable version of sgRSEA] (http://)
[Download the R package of sgRSEA] (https://cran.r-project.org/web/packages/sgRSEA/index.html)

For detailed [manual](wiki page)

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
'''python
python setup.py install
```

##Usage






##File formats
###Design file
When you are using the suite for your experiment, you need to prepare a design file.(For individual functionalities, you may not have to.)
Design file has to include all the essential columns (Please use exactly the same column names, order does not matter):

* **filepath**: The absolute path to the fastq file
* **lib**: The absolute path the the library file. If there is no sublib, this column should has the same value across all rows
* **sublib**: The label for sublib. Eg: GeCKO_libA, GeCKO_libB. Use this when you sequence sublib separately
* **label**: The user friendly label for each sample. This will be used as output prefix for each sample. Please DON'T use " ", "-" in the label 
* **sgstart**: The first nucleotide of sgRNA. 1-index
* **sgstop**: The last nucleotide of sgRNA. 1-index
* **trim3**: Sequence pattern of the 3' adaptor. Usually 5~7nt. If provided, the program will look for perfect match of this pattern in fastq sequence. The last match and all nucleotides after that will be trimmed. If you don't need this, put "NA" in the design file


 


















https://cran.r-project.org/web/packages/sgRSEA/index.html