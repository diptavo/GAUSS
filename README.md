# About

GAUSS: **G**ene-set **A**ssociation analysis **U**sing **S**parse **S**ignals

GAUSS provides a powerful and computationally effective tool to perform gene-set (pathway) association analysis. The methods tests for any self-contained association between a phenotype and a gene-set and produces a p-value for the association using results from gene-based tests (e.g. SKAT, SKAT-O, SKAT-Common-Rare, Burden test, prediXcan etc.).In addition, if there is a significant association, GAUSS identifies the genes within the gene-set that might be driving the association.

# Citation 

If you GAUSS to analyze gene-set associations, please consider citing our [biorXiv paper](https://www.biorxiv.org/content/10.1101/799791v2). 

If you have any questions/issues about the package, please email me at **diptavo21@jhu.edu**

# Input files

GAUSS requires two input files

- Summary gene-based p-value file in flat-text format, **without** headers. The file should contain at least two columns: Gene-name and the corresponding p-values. An example has been provided: *example_gene_pval.txt*. The columns do not need to be in a particular order and there can be additional columns as well. Please use **Gene-Symbols** for the genes.

- GMT file in flat text format, **with** headers. This file should be in the format as specified in *example_gmt.txt*. The first column contains the name of the gene-set, second column contains information on it (possibly URLs) and the thrid column contains the list of genes in the gene-set comma separated.

# Prerequisites and download

GAUSS was built using R (v 3.6.0). In general any version of R `>= 3.5` should work. 
The GAUSS codes depend on several existing R packages: `optparse`, `data.table`, `gPdtest`, `POT`, `MASS`. Please download and install them prior to running GAUSS.

The repository can be cloned and installed as:

```
 git clone https://github.com/diptavo/GAUSS.git
 cd GAUSS/
 R CMD INSTALL GAUSS_1.0.tar.gz 
```

# Usage

Several example files are provided with the package for the user to verify the formats and for toy-examples.
- GMT files: `example_gmt.txt` and `example_gmt2.txt`
- gene-based summary p-value files: `example_gene_pval.txt`

Once installed, GAUSS can be run with following commands (assuming the path to GAUSS repository is `~/GAUSS/`)

```R

library(GAUSS)

GAUSS_All(summary_file = "~/GAUSS/example_gene_pval.txt", gene_name = 1, pv_name = 2, output_file = "example_out", gmt = "~/GAUSS/example_gmt.txt", ags = "def",verbose = TRUE,parallel = FALSE)
### This will produce two files: example_out.log and example_out.out in about 0.7 minutes. 

GAUSS_All(summary_file = "~/GAUSS/example_gene_pval.txt", gene_name = 1, pv_name = 2, output_file = "example_out", gmt = "~/GAUSS/example_gmt2.txt", ags = "def",verbose = TRUE,parallel = FALSE)
### This will produce two files: example_out.log and example_out.out in about 6 minutes.  
```
The output files will contain the following information:
- `example_out.out`: Flat text file containing gene-set, GAUSS p-value and selected CS genes in an R readable format using `read.table()`.
- `example_out.log`: Log file containing run-time information for each gene-set and the overall run-time.

To run GAUSS on a subset of gene-sets present in the GMT file, use the following command:

```R
GAUSS_All(summary_file = "~/GAUSS/example_gene_pval.txt", gene_name = 1, pv_name = 2, output_file = "example_out", gmt = "~/GAUSS/example_gmt2.txt", ags = "def",verbose = TRUE,parallel = TRUE,start = 21,stop = 50)

### This will run GAUSS on the 21st gene-set (GO_FOREBRAIN_NEURON_DEVELOPMENT) through 50th gene-set (GO_ACYLGLYCEROL_HOMEOSTASIS)
```

## Options for `GAUSS_All(.)`

- `summary_file`: Summary file containing the gene-names and corresponding p-values

- `gene_name`: Column number for gene names in summary file

- `pv_name`: Column number for p-value in summary file

- `output_file`: prefix for output file names

- `gmt`: GMT file containing the list of gene-set, one gene-set per line

- `verbose`: Print extra output; default = TRUE

- `ags`: settings for control arguments for running GAUSS; options are "def" or "prec"; "def" should be used for a quick initial scan and the significant associations can be followed up using "prec".

- `parallel`: logical indicating whether parallel jobs are to be created; default = FALSE

- `jobs`: number of parallel jobs to be created if parallel = TRUE

- `jobfile`: name of the output file with parallel jobs for GAUSS if parallel = TRUE

- `start`: starting point of the GMT file if parallel = TRUE

- `stop`: stopping point of the GMT file if parallel = TRUE




# Run parallel jobs

For a GMT file containing many gene-set definitions, it is much easier to divide it into multiple jobs as:

```shell

Rscript ~/GAUSS/utils/run_GAUSS_All.R --summary ~/GAUSS/example_gene_pval.txt --gmtFile ~/GAUSS/example_gmt.txt --out ex1 --pvalue 2 --geneName 1 --parallel TRUE --jobs 10 -f ~/AJ.txt

```
This will create 10 `Rscript` jobs in the file `~/AJ.txt`. 

```shell
head -2 ~/AJ.txt

Rscript run_GAUSS_All.R --summary /Users/diptavo/GAUSS/example_gene_pval.txt --geneName 1 --pvalue 2 --out ex1_1 --gmtFile /Users/diptavo/GAUSS/example_gmt.txt --verbose TRUE --ags def --parallel FALSE --start 1 --stop 20
Rscript run_GAUSS_All.R --summary /Users/diptavo/GAUSS/example_gene_pval.txt --geneName 1 --pvalue 2 --out ex1_2 --gmtFile /Users/diptavo/GAUSS/example_gmt.txt --verbose TRUE --ags def --parallel FALSE --start 21 --stop 40

```

You can run each line of this file appropriately in batch mode, HPC etc.  

## Options for `run_GAUSS_All.R`

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

# Update Log

- 07/06/2020: GAUSS published as a R-package along with detailed documentations and utility function.
