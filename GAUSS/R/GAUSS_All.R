GAUSS_All <- function(summary_file, gene_name, pv_name,output_file,gmt,verbose = TRUE,ags = "def",parallel = FALSE,jobs,jobfile,start,stop,is.appx = TRUE){
  
  LOGF<- paste0(output_file,".log")
  OUTF <-paste0(output_file,".out") 
  ll <- file(LOGF,open = "wt")
  sink(ll,type = "message")
  
  cat("Starting GAUSS: \n\n",file = LOGF,append = T)
  cat("options parsed \n\n",file = LOGF,append = T)
  if(is.appx){
    cat(paste0("GPD approximations for small p-values will be used ... \n"),file = LOGF, append = T)
  }else{
    cat(paste0("No approximations invoked. Only resampling p-values of limited precision will be obtained... \n"),file = LOGF,append = T)
  }
  print("loading ref data...")
  cat("loading ref data ...\n\n",file = LOGF,append = T)
  print(head(colnames(pv.null.wt1)))
  
  if(ags == "def"){
    ags = c(1e+05,250,100,-10,500,0.05)
  }else if(ags == "prec"){
    ags = c(1e+06,250,100,-10,500,0.05)
  }else{
    ags <- read.table(ags)[,1]
  }
  
  s1 <- read.table(gmt,header = T); summary_file = read.table(summary_file);
  if(parallel){
    start <- start; stop <- stop; if(is.na(start) || is.na(stop)){stop("Need a start and stop index for GMT file when running in parallel")}
  }else{
    if(is.na(start) || is.na(stop)){
      start <- 1; stop <- dim(s1)[1]
    }
  }
  
  s1 <- s1[c(start:stop),]
  print(paste0("reading summary and gmt files... "))
  cat("reading summary and gmt files \n\n",file = LOGF,append = T)
  pv <- NULL;
  cs <- list();
  set.seed(3721)
  
  print(paste0("logging in ...",LOGF))
  print(paste0("output in ... ",OUTF))
  
  
  cat(paste0("logging in ... ",LOGF,"\n\n"),file = LOGF,append = T)
  cat(paste0("output in ... ",OUTF,"\n\n"),file = LOGF,append = T)
  
  
  aaa <- system.time(for(i in 1:dim(s1)[1]){
    
    if(verbose){
      cat(paste0("Running gene-set ",i,": ",s1[i,1]),file = LOGF,append = T)
      cat(paste0("\n"), file = LOGF,append = T)
    }
    
    rt1 <- system.time(t1 <- gauss(gene.set= as.character(s1[i,3]),inf = summary_file,ggg = gene_name,cols1 = pv_name,max.prec = ags[1], ct=ags[2],ct.min = ags[3],stp = ags[4],tl = ags[5],thr = ags[6],verbose = verbose,LOGF = LOGF,is.appx = is.appx))
    pv <- c(pv,t1$pvalue); cs[[i]] <- t1$CS;
    cat(paste0("Analysis done in ",round(as.numeric(rt1[3]),4)," seconds \n\n"),file = LOGF, append = T)
    write.table(data.frame(s1[1:i,1],pv[1:i],I(unlist(lapply(cs,paste,collapse=",")))),OUTF,col.names = c("GeneSet","pvalue","CS"),row.names = F,quote = F);
  })
  
  cat(paste0("\n\n"), file = LOGF,append = T)
  cat(paste0("Analyzed ",i," gene-sets in ",round(as.numeric(aaa[3]/60),5)," minutes \n\n"), file = LOGF, append = T)
  sink(type = "message")
  close(ll)
  
}