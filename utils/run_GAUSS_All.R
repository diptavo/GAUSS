# args = commandArgs(trailingOnly=TRUE)
library(optparse)
library(data.table)
library(gPdtest)
library(POT)
library(MASS)


option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-s", "--summary"), type="character", 
              help="Summary file containing the gene-names and corresponding p-values [default %default]"),
  make_option(c("-n", "--geneName"), type = "integer", default= 1, 
              help = "Column number for gene names in summary file [default %default]"),
  make_option(c("-g", "--gmtFile"), type = "character", 
              help = "GMT file containing the list of gene-set, one gene-set per line"),
  make_option(c("-p", "--pvalue"), default=2, 
              help="Column number for p-value in summary file [default %default]"),
  make_option(c("-a", "--ags"),type = "character", default = "def",
              help="control arguments for running GAUSS [default %default]"),
  make_option(c("-o","--out"),type = "character",default="out",
              help="prefix for output file names [default %default]"),
  make_option(c("-d","--path"),type = "character",default="~/GAUSS",
              help="path to files for GAUSS [default %default]"),
  make_option(c("-r","--parallel"),default=FALSE,
              help="logical indicating whether parallel jobs are to be created [default %default]"),
  make_option(c("-j","--jobs"),type = "integer",default=1,
              help="number of parallel jobs to be created [default %default]"),
  make_option(c("-f","--parallelFile"),type = "character",default = "jobs.txt",
              help="name of the output file with parallel jobs for GAUSS [default %default]"),
  make_option(c("-t", "--start"),type = "double",default = NA,
              help="the starting point in GMT file if a parallel job is invoked [default %default]"),
  make_option(c("-z", "--stop"),type = "double",default = NA,
              help="the stopping point in GMT file if a parallel job is invoked [default %default]")
  
)
opt <- parse_args(OptionParser(option_list=option_list))

summary_file = as.character(opt$summary)
gene_name = as.integer(opt$geneName)
pv_name = as.integer(opt$pvalue)
output_file = as.character(opt$out)
gmt <- as.character(opt$gmtFile)
verbose <- as.character(opt$verbose)
ags <- as.character(opt$ags)
is.parallel <- opt$parallel
jobs <- as.integer(opt$jobs)
jobfile <- as.character(opt$parallelFile)
start <- as.numeric(opt$start)
stop <- as.numeric(opt$stop)

cat(paste0("Starting GAUSS...\n\n"))

library(GAUSS)

if(is.parallel){
sink(jobfile)
cat(paste0("trying to open file..."))
sink()
jobfile <- normalizePath(jobfile)
}


if(is.parallel){
  s1 <- read.table(gmt,header = T); npath <- nrow(s1); ind <- ceiling(npath/jobs); indx <- seq(0,npath,ind)
  if(npath%%jobs != 0){
    sink(jobfile)
    for(i in 1:(jobs-1)){
      OUTF <- paste0(output_file,"_",i)
      cat(paste0("Rscript run_GAUSS_All.R --summary ",summary_file," --geneName ",gene_name," --pvalue ",pv_name," --out ",OUTF," --gmtFile ",gmt," --verbose ",verbose," --ags ",ags," --parallel FALSE --start ",indx[i]+1," --stop ",indx[i+1]))
      cat(paste0("\n"));
    }
    cat(paste0("Rscript run_GAUSS_All.R --summary ",summary_file," --geneName ",gene_name," --pvalue ",pv_name," --out ",OUTF," --gmtFile ",gmt," --verbose ",verbose," --ags ",ags," --parallel FALSE --start ",indx[i+1]+1," --stop ",npath))  
    sink()
  }else{
    sink(jobfile)
    for(i in 1:(jobs)){
      OUTF <- paste0(output_file,"_",i)
      cat(paste0("Rscript run_GAUSS_All.R --summary ",summary_file," --geneName ",gene_name," --pvalue ",pv_name," --out ",OUTF," --gmtFile ",gmt," --verbose ",verbose," --ags ",ags," --parallel FALSE --start ",indx[i]+1," --stop ",indx[i+1]))
      cat(paste0("\n"));
    }
    sink()
  }
  cat(paste0(npath," pathway definitions in GMT files.. \n\n")); jlen <- system(paste0("wc -l ",normalizePath(jobfile)),intern = T)
  cat(paste0(jlen," parallel jobs are created in ",jobfile," .... \n\n"))
  cat(paste0("See ",jobfile,"... You can appropriately run each line as a batch job \n\n"))
}

if(!is.parallel){
  GAUSS_All(summary_file = summary_file,
            gene_name = gene_name,
            pv_name = pv_name,
            output_file = output_file,
            gmt = gmt,
            verbose = verbose,
            ags = ags,
            parallel = FALSE,
            start = start,
            stop = stop)}
