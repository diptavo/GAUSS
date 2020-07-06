Calc.PV <-
  function(stat1,cor.z,max.prec = 1e+05,init = 1000,step = 100,oc = 2,ct=250,ct.min = 100,stp = -10,tl = 500,thr = 0.05,verbose = T,LOGF = "out.log",is.appx = TRUE){
    
    if(verbose)
      cat("Running resampling \n",file = LOGF,append = T)
    precs <- init; gpd.used = FALSE;r <- Resample.PV(stat1,cor.z,prec = precs)
    pval <- r$pval;max.stat <- r$max.stat
    while(pval < (2/init) && precs < max.prec){
      precs <- precs*step
      r <- Resample.PV(stat1,cor.z,prec = precs)
      pval <- r$pval; max.stat <- r$max.stat
    }
    if(pval < oc/length(max.stat)){
      if(is.appx){
        pval <- gpd.pv(max.stat,tst = stat1,ct=ct,ct.min = ct.min,stp = stp,tl = tl,thr = thr,LOGF = LOGF);
        gpd.used = TRUE
      }
    }
    
    return(list("pval" = pval,"gpd.used" = gpd.used))
  }
