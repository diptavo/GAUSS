# About

GAUSS: **G**ene-set **A**ssociation analysis **U**sing **S**parse **S**ignals

GAUSS provides a powerful and computationally effective tool to perform gene-set (pathway) association analysis. The methods tests for any self-contained association between a phenotype and a gene-set and produces a p-value for the association using results from gene-based tests (e.g. SKAT, SKAT-O, SKAT-Common-Rare, Burden test, prediXcan etc.).In addition, if there is a significant association, GAUSS identifies the genes within the gene-set that might be driving the association.

# Citation 

If you use GAUSS, please cite our article now publishd in *AJHG*: [A powerful subset-based method identifies gene set associations and improves interpretation in UK Biobank](https://www.sciencedirect.com/science/article/abs/pii/S0002929721000586)

You can find a preprint version of the paper on [biorXiv](https://www.biorxiv.org/content/10.1101/799791v2). 

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

If you face trouble in installing the tar zipped file in the above procedure, please download the file from this [link](https://www.dropbox.com/s/iq6w139fag3ozvw/GAUSS_1.0.tar.gz?dl=1)


# Usage

Several example files are provided with the package for the user to verify the formats and for toy-examples.
- GMT files: `example_gmt.txt` and `example_gmt2.txt`
- gene-based summary p-value files: `example_gene_pval.txt`

Once installed, GAUSS can be run with following commands (assuming the path to GAUSS repository is `~/GAUSS/`)

```R

library(GAUSS)

GAUSS_All(summary_file = "~/GAUSS/example_gene_pval.txt", gene_name = 1, pv_name = 2, output_file = "example_out", gmt = "~/GAUSS/example_gmt.txt", ags = "def",verbose = TRUE,parallel = FALSE)
closeAllConnections()
### This will produce two files: example_out.log and example_out.out in about 0.7 minutes. 

GAUSS_All(summary_file = "~/GAUSS/example_gene_pval.txt", gene_name = 1, pv_name = 2, output_file = "example_out", gmt = "~/GAUSS/example_gmt2.txt", ags = "def",verbose = TRUE,parallel = FALSE)
closeAllConnections()
### This will produce two files: example_out.log and example_out.out in about 6 minutes.  
```
The output files will contain the following information:
- `example_out.out`: Flat text file containing gene-set, GAUSS p-value and selected CS genes in an R readable format using `read.table()`.
- `example_out.log`: Log file containing run-time information for each gene-set and the overall run-time.

To run GAUSS on a subset of gene-sets present in the GMT file, use the following command:

```R
GAUSS_All(summary_file = "~/GAUSS/example_gene_pval.txt", gene_name = 1, pv_name = 2, output_file = "example_out", gmt = "~/GAUSS/example_gmt2.txt", ags = "def",verbose = TRUE,parallel = TRUE,start = 21,stop = 50)
closeAllConnections()

### This will run GAUSS on the 21st gene-set (GO_FOREBRAIN_NEURON_DEVELOPMENT) through 50th gene-set (GO_ACYLGLYCEROL_HOMEOSTASIS)
```

## Options for `GAUSS_All(.)`

- `summary_file`: Summary file containing the gene-names and corresponding p-values

- `gene_name`: Column number for gene names in summary file

- `pv_name`: Column number for p-value in summary file

- `output_file`: prefix for output file names

- `gmt`: GMT file containing the list of gene-set, one gene-set per line

- `verbose`: Print extra output; default = TRUE

- `method`: Method/test used to generate the gene-based p-values. The reference Vh correlation will be loaded according to this. Currently the supported options are "SKAT-CR" (for SKAT Common-Rare) or "TWAS.X" (for TWAS on tissue X. e.g: "TWAS.Whole_Blood"). For full list of available reference Vh, see `data/` directory. 

- `ags`: settings for control arguments for running GAUSS; options are "def" or "prec"; "def" should be used for a quick initial scan and the significant associations can be followed up using "prec".

- `parallel`: logical indicating whether parallel jobs are to be created; default = FALSE

- `jobs`: number of parallel jobs to be created if parallel = TRUE

- `jobfile`: name of the output file with parallel jobs for GAUSS if parallel = TRUE

- `start`: starting point of the GMT file if parallel = TRUE

- `stop`: stopping point of the GMT file if parallel = TRUE

- `is.appx`: GPD approximation invoked. default = TRUE

- `gpd.est.res`: Number of resampling iterations to estimate GPD parameters. default = 250

- `pv.null.wt1`: User input reference dataset for calculating Vh. Please refer to the manuscript and the next section for more details.

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

- `--verbose` : Outputs some run-time messages in `.log` file; default is `TRUE`
- `--ags` : A path to a file that contains some control options for GAUSS. Please don't change this. This will be updated in a later version.

For full list of options please run

```
Rscript ~/GAUSS/utils/run_GAUSS_All.R --help

```
# Results

We performed the association analysis of 1,403 binary phenotypes from UK-Biobank with `C2` (Curated pathways) and `C5` (GO pathways) from [MSigDB v6.2](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/6.2/) using GAUSS. The results can be visualized using a [PheWeb-like visual server](http://ukb-pathway.leelabsg.org/). 

## Creating own reference data (Vh) for customized gene-based tests

You can create your own reference data for conducting GAUSS test if you are not using any of the standard tests included in GAUSS. To create a compatible Vh dataset follow these steps:

- You need the individual level genotype files and an annotation file. In absence of individual level genotypes, an acceptable set of GWAS variants can be obtained from the 1000 Genomes genotypes obtained via the [FUSION or LDSCore packages](http://gusevlab.org/projects/fusion/#installation) `wget https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2`

- The annotation file would contain information on the assignment of variants to genes. The format needs to be compatible with the gene-based test software you are going to use. 

- Independently generate a null phenotype from `N(0,1)` without any variant effects.

- Perform gene-based test using the test/software of your choice, annotation files and individual level genotype data and simulated phenotype data.

- Record the p-values of the genes and convert them to z-values.

- Repeat this step multiple times (100 or 500 or 1000) and obtain a matrix of z-values with dimensions `RxG` where `R` is the number of repeatations (= 100/500/1000 or others as specified) and `G` is the number of genes. Store this matrix and use this as input with option `pv.null.wt1` in `GAUSS_All(.)`

An example workflow to create TWAS-FUSION reference Vh for a particular tissue:

- Independently generate a null phenotype from `N(0,1)` without any variant effects.

- Run GWAS with simulated phenotype and variants from 1000 Genomes obtained at the (FUSION package)[https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2] using GWAS pipelines like PLINK, GCTA/fastGWA or EPACTS

- For a given tissue, run TWAS analysis using the summary statistics from the above GWAS and appropriate tissue-specific weights following the (instructions)[http://gusevlab.org/projects/fusion/#typical-analysis-and-output]

- Repeat this multiple (`R`) times and store the resultant null z-values in a matrix format as described above. 


# Update Log

- 06/13/2021: Bugs in options and parallel job creation fixed
- 01/04/2021: README updated
- 12/16/2020: TWAS-FUSION reference panel added
- 07/06/2020: GAUSS published as a R-package along with detailed documentations and utility function.
