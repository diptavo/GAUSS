Resample.PV.C <-
function(stat1,cor.z,condMean,prec = 10000){
  
  
  max.stat <- NULL
  b1<-Re(mat.sqrt(cor.z))
  ran1<-matrix(rnorm(dim(cor.z)[1]*prec),ncol=prec)
  ran2<- t(b1)%*% ran1
  system.time(for(id in 1:prec){
    a1<-ran2[,id] + condMean
    max.stat <- c(max.stat,Cal_S_EqualWeight(a1)$maxS)
  })
  
  pval <- sum(max.stat > stat1)/length(max.stat)
  return(pval)
}
