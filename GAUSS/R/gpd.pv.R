gpd.pv <-
function(max.stat,tst,ct=250,ct.min = 100,stp = -10,tl = 500,thr = 0.05,LOGF = "out.log"){
  
  cat("approximation using GPD \n\n",file = LOGF,append = T)
  cutoff <- ct
  sq <- seq(ct,ct.min,by = stp)
  for(s in sq){
    cutoff <- s
    ms1 <- max.stat[order(max.stat,decreasing = T)]; t <- ms1[(cutoff+1)]
    ms1 <- ms1[1:cutoff] - t
    test.obj1 <- gpd.test(ms1,J=tl)
    if(test.obj1$boot.test$p.value > thr){
      break;
    }
  }
  
  ss = tst
  f2 = fitgpd(max.stat,threshold=t,est = "mgf",stat = "AD2R")
  pv3 <- pgpd((ss-t),shape = f2$fitted.values[2],scale = f2$fitted.values[1],lower.tail = F)
  pv.gpd <- as.numeric(pv3*s/length(max.stat))
  if(pv.gpd < 1e-200){
    pv.gpd <- pexp(ss-t,rate = 1/mean(ms1),lower.tail = F)*s/length(max.stat)
  }
  return(pv.gpd)
}
