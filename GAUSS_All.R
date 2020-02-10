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
              help="path to files for GAUSS [default %default]")
)
opt <- parse_args(OptionParser(option_list=option_list))

summary_file = as.character(opt$summary)
gene_name = as.integer(opt$geneName)
pv_name = as.integer(opt$pvalue)
output_file = as.character(opt$out)
gmt <- as.character(opt$gmtFile)
verbose <- as.character(opt$verbose)
ags <- as.character(opt$ags)
path <- as.character(opt$path)

LOGF<- paste0(output_file,".log")
OUTF <-paste0(output_file,".out") 
ll <- file(LOGF,open = "wt")
sink(ll,type = "message")

cat("options parsed \n\n",file = LOGF,append = F)

cat("Starting GAUSS: \n\n",file = LOGF,append = T)
source(paste0(path,"/GAUSS.R")); source(paste0(path,"/GPD.R")); source(paste0(path,"/Helper.R"));
print("loading ref data...")
cat("loading ref data ...\n\n",file = LOGF,append = T)
load("~/GAUSS/Null.RData")

if(ags == "def"){
  ags = c(250,100,-10,500,0.05)
}else{
  ags <- read.table(ags)[,1]
}

# cl.sumstat <- summary(file(summary_file))$class;

# if(cl.sumstat == "gzfile"){
#   print("File in gzip format.. \n")
#   print("Reading gene-based test files in .gz format... \n")
#   t1 <- fread(paste("zcat ",summary_file))
# }else{   
#   print("Reading gene-based test files in plain text format... \n")
#   t1 <- fread(paste(summary_file))
#   print("File not in gzipped format")
# }        

s1 <- read.table(gmt,header = T); summary_file = read.table(summary_file);
cat("reading summary and gmt files \n\n",file = LOGF,append = T)
pv <- NULL;
cs <- list();
print(paste0("logging in ...",LOGF))
cat(paste0("logging in ...",LOGF,"\n\n"),file = LOGF,append = T)

aaa <- system.time(for(i in 1:dim(s1)[1]){
  
  if(verbose){
    cat(paste0("Running gene-set ",i),file = LOGF,append = T)
    cat(paste0("\n"), file = LOGF,append = T)
  }
  
  t1 <- gauss(gene.set= as.character(s1[i,3]),inf = summary_file,ggg = gene_name,cols1 = pv_name,ct=ags[1],ct.min = ags[2],stp = ags[3],tl = ags[4],thr = ags[5],verbose = verbose,LOGF = LOGF)
  pv <- c(pv,t1$pvalue); cs[[i]] <- t1$CS;
  write.table(data.frame(s1[1:i,1],pv[1:i],I(unlist(lapply(cs,paste,collapse=",")))),OUTF,col.names = c("GeneSet","pvalue","CS"),row.names = F,quote = F);
})

cat(paste0("\n\n"), file = LOGF,append = T)
cat(paste0("Analyzed ",i," gene-sets in ",round(as.numeric(aaa[3]/60),5)," minutes \n\n"), file = LOGF, append = T)
sink(type = "message")
close(ll)
