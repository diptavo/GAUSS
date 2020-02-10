# About

GAUSS: **G**ene-set **A**ssociation analysis **U**sing **S**parse **S**ignals

GAUSS provides a powerful and computationally effective tool to perform gene-set (pathway) association analysis. The methods tests for any self-contained association between a phenotype and a gene-set and produces a p-value for the association using results from gene-based tests (e.g. SKAT, SKAT-O, SKAT-Common-Rare, Burden test, prediXcan etc.).In addition, if there is a significant association, GAUSS identifies the genes within the gene-set that might be driving the association.

# Citation 

If you GAUSS to analyze gene-set associations, please consider citing our [biorXiv paper](https://www.biorxiv.org/content/10.1101/799791v1). 

If you have any questions/issues about the package, please email me at **diptavo21@jhu.edu**

# Input files

GAUSS requires two input files

- Gene-based p-value file in flat-text format, **without** headers. The file should contain at least two columns: Gene-name and the corresponding p-values. An example has been provided: *example_gene_pval.txt*

- GMT file in flat text format, **with** headers. This file should be in the format as specified in *example_gmt.txt*. The first column contains the name of the gene-set, second column contains information on it (possibly URLs) and the thrid column contains the list of genes in the gene-set comma separated.

# Prerequisites and download

The GAUSS codes depend on several existing R packages: *optparse, data.table, gPdtest, POT, MASS*. Please download and install them prior to running GAUSS.

The repository can be cloned as:

'git clone https://github.com/diptavo/GAUSS.git'

# Usage

Once downloaded, GAUSS can be run with following commands (assuming the path to GAUSS is `~/GAUSS/`)

`Rscript --vanilla ~/GAUSS/GAUSS_All.R --summary ~/GAUSS/example_gene_pval.txt --gmtFile ~/GAUSS/example_gmt.txt --out ex1 --pvalue 2 --geneName 1 &`

