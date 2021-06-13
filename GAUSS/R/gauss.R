gauss <-
  function(gene.set,inf,ggg,cols1,max.prec = 1e+05, ct=250,ct.min = 100,stp = -10,tl = 500,thr = 0.05,verbose = T,LOGF = "out.log",is.appx = T){
    
    inf = inf[match(colnames(pv.null.wt1),inf[,ggg]),]; inf <- inf[!is.na(inf[,cols1]),]  
    gene.set <- gsub(c(" "), "",gene.set); gene.set <- unique(gene.set)
    gene.set <- strsplit(as.character(gene.set),split = ",")[[1]]
    length.path <- length(gene.set)
    if(verbose)
      cat("Preprocessing \n",file = LOGF,append = T)
    
    
    common <- intersect(gene.set,colnames(pv.null.wt1)); excl <- length.path - length(common);
    if(length(common) > 0){
      gs1 <- gene.set[match(common,gene.set)]
      cl1 <- colnames(pv.null.wt1)[match(common,colnames(pv.null.wt1))]
      
      z.mat <- as.matrix(pv.null.wt1[,colnames(pv.null.wt1) %in% cl1])
      
      if(length(z.mat) > 0){
        cor.z = cor(as.matrix(z.mat),use = "pairwise.complete.obs");
        
        pv.mat <- inf[match(common,colnames(pv.null.wt1)),]
        pv1 <- pv.mat[,cols1]
        if(sum(is.na(pv1)) == length(pv1)){
          return(list("pvalue" = NA, "CS" = NA, "Test.stat" = NA,"gpd.used" = NA,"excl" = excl, "length" = common))
        }else{
          p.noinf = pv1 != 1
          
          if(sum(p.noinf,na.rm = T) > 0){
            pv2 = pv1[p.noinf]
            cor.z = as.matrix(cor.z[p.noinf,p.noinf])
            
            pv.mat2 <- pv.mat[p.noinf,]
            diag(cor.z) = 1; cor.z[is.na(cor.z)] = 0; cor.z[is.nan(cor.z)] = 0;
            
            if(sum(p.noinf,na.rm = T) > 0){
              test.obj <- Cal_S_EqualWeight(qnorm(pv2,lower.tail = F))
              stat.wt <- test.obj$maxS
              gene.idx.wt <- pv.mat2[test.obj$idx,1]
              p1 <- Calc.PV(stat.wt,cor.z,max.prec = max.prec, ct=ct,ct.min = ct.min,stp = stp,tl = tl,thr = thr,verbose = verbose,LOGF = LOGF,is.appx = is.appx)
              pval.wt <- p1$pval; gpd.used = as.numeric(p1$gpd.used)
            }
            
            gauss.obj <- list("pvalue" = pval.wt, "CS" = gene.idx.wt, "Test.stat" = stat.wt,"gpd.used" = gpd.used, "excl" = excl, "length" = length.path);
            return(gauss.obj)
          }else{
            return(list("pvalue" = NA, "CS" = NA, "Test.stat" = NA,"gpd.used" = NA,"excl" = excl, "length" = common))
          }
        }
      }else{
        return(list("pvalue" = NA, "CS" = NA, "Test.stat" = NA,"gpd.used" = NA,"excl" = excl, "length" = common))
      }
      
    }else{
      return(list("pvalue" = NA, "CS" = NA, "Test.stat" = NA,"gpd.used" = NA,"excl" = excl, "length" = common))
    }
  }
