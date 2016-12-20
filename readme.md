#sgRSEA
###Enrichment Analysis of CRISPR/Cas9 Knockout Screen Data




**Version** : v0.0.5

**Authors**

Beibei Chen (beibei.chen@utsouthwestern.edu)
Jungsik Noh (junsik.noh@utsouthwestern.edu)

**Maintainer**: Beibei Chen

For detailed [manual](https://github.com/bchen4/sgrsea/wiki/Manual)


##Prerequisites
* Python 2.7
* [Pandas](http://pandas.pydata.org/)
* [Numpy](http://www.numpy.org/)
* [statsmodels](http://statsmodels.sourceforge.net/) 
* [scipy](https://www.scipy.org/)
* [biopython](http://biopython.org/)

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

##Quick start

sgRSEA can be run on specific steps as well as fastq-to-result using a single command.

```python
#Run the suite from fastq to final results
sgrsea run -d design_file -o my_experiment

#Get the count for single fastq file
sgrsea count -i my_fastq -o output -l sgRNA_lib 

#Get the count matrix of multiple fastq file and generate a matrix
sgrsea count -d design.txt -o output_prefix

#Normalize the count matrix
sgrsea normalization -i count_matrix --ormalize-method total -o outputfile

#Normalize the count matrix. There are sub-libs sequenced separately
sgrsea normalization -i count_matrix --ormalize-method total -o outputfile --split-lib

#Reformat the matrix into sgRSEA 4-column matrix. 
#-t and -c values MUST match group content of the design file
sgrsea reformat -i normalized_matrix -d design.txt -t Heat,Cold,Dry -c Ctrl,Ctrl,Ctrl -o output_prefix --collapse-replicates auto

#Stattest on normalized and formatted matrix
sgrsea stattest -i matrix -o output


```

##License
See the [LICENSE](https://github.com/bchen4/sgrsea/blob/master/LICENSE.txt) file for license rights and limitations MIT.


##Also available in R

https://cran.r-project.org/web/packages/sgRSEA/index.html
