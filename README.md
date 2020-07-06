# About

GAUSS: **G**ene-set **A**ssociation analysis **U**sing **S**parse **S**ignals

GAUSS provides a powerful and computationally effective tool to perform gene-set (pathway) association analysis. The methods tests for any self-contained association between a phenotype and a gene-set and produces a p-value for the association using results from gene-based tests (e.g. SKAT, SKAT-O, SKAT-Common-Rare, Burden test, prediXcan etc.).In addition, if there is a significant association, GAUSS identifies the genes within the gene-set that might be driving the association.

# Citation 

If you GAUSS to analyze gene-set associations, please consider citing our [biorXiv paper](https://www.biorxiv.org/content/10.1101/799791v1). 

If you have any questions/issues about the package, please email me at **diptavo21@jhu.edu**

# Input files

GAUSS requires two input files

- Gene-based p-value file in flat-text format, **without** headers. The file should contain at least two columns: Gene-name and the corresponding p-values. An example has been provided: *example_gene_pval.txt*. The columns do not need to be in a particular order and there can be additional columns as well. Please use **Gene-Symbols** for the genes.

- GMT file in flat text format, **with** headers. This file should be in the format as specified in *example_gmt.txt*. The first column contains the name of the gene-set, second column contains information on it (possibly URLs) and the thrid column contains the list of genes in the gene-set comma separated.

# Prerequisites and download

GAUSS was built using R (v 3.6.0).The GAUSS codes depend on several existing R packages: `optparse`, `data.table`, `gPdtest`, `POT`, `MASS`. Please download and install them prior to running GAUSS. Upon installing the required packages no further installation is necessary.

The repository can be cloned as:

```
 git clone https://github.com/diptavo/GAUSS.git
 cd GAUSS/
 R CMD INSTALL GAUSS_1.0.tar.gz 
```

# Usage

Several example files are provided with the package for the user to verify the formats and for toy-examples.
Once installed, GAUSS can be run with following commands (assuming the path to GAUSS repository is `~/GAUSS/`)

```R

library(GAUSS)

GAUSS_All(summary_file = "~/GAUSS/example_gene_pval.txt", gene_name = 1, pv_name = 2, output_file = "example_out", gmt = "~/GAUSS/example_gmt.txt", ags = "def",verbose = TRUE,parallel = FALSE)

GAUSS_All(summary_file = "~/GAUSS/example_gene_pval.txt", gene_name = 1, pv_name = 2, output_file = "example_out", gmt = "~/GAUSS/example_gmt2.txt", ags = "def",verbose = TRUE,parallel = FALSE)

### This will produce two files: example_out.log and example_out.out in about 0.7 minutes.  
```
The output files will contain the following information:
`example_out.out`: Flat text file containing gene-set, GAUSS p-value and selected CS genes in an R readable format using `read.table()`.
`example_out.log`: Log file containing run-time information for each gene-set and the overall run-time.

# Options

- `--summary` : Gene-based p-value file in flat-text format **without** header.
- `--geneName` : the column number in the above file that contains the names of the genes; defualt value is `1`. 
- `--pvalue` : the column number in the above file that contains the p-values for the genes; defualt value is `2`. 
- `--out`: the prefix of the output file; default is `out`. Two files will be generated: `.log` contains some information on running and `.out` contains the gene-set, p-value and the corresponding CS genes.
- `--gmtFile` : GMT file containing the list of gene-set, one gene-set per line in the format mentioned above, **with** a header line.

There are a few other options that can be used:

- `--path` : The path to the GAUSS directory. By default it is set to `~/GAUSS/`
- `--verbose` : Outputs some run-time messages in `.log` file; default is `TRUE`
- `--ags` : A path to a file that contains some control options for GAUSS. Please don't change this. This will be updated in a later version.

# Results

We performed the association analysis of 1,403 binary phenotypes from UK-Biobank with `C2` (Curated pathways) and `C5` (GO pathways) from [MSigDB v6.2](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/6.2/) using GAUSS. The results can be visualized using a [PheWeb-like visual server](http://ukb-pathway.leelabsg.org/). 

# Miscellaneous 

This is an initial developement version. This will be updated along with the manuscript. 
