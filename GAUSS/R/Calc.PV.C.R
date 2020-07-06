Calc.PV.C <-
function(stat1,cor.z,condMean,max.prec = 1e+05,init = 1000,step = 100,LOGF = "out.log"){
  
  cat("Running resampling",file = LOGF,append = T)
  precs <- init
  pval <- Resample.PV.C(stat1,cor.z,condMean = condMean,prec = precs)
  while(pval < (2/init) && precs < max.prec){
    precs <- precs*step
    pval <- Resample.PV(stat1,cor.z,prec = precs)
  }
  return(pval)
}
