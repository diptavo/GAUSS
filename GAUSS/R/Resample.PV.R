Resample.PV <-
function(stat1,cor.z,prec = 100000){
  
  max.stat <- NULL
  b1<-mat.sqrt(cor.z)
  ran1<-matrix(rnorm(dim(cor.z)[1]*prec),ncol=prec)
  ran2<- t(b1)%*% ran1
  max.stat <- apply(ran2,2,CMS)
  # system.time(for(id in 1:prec){
  #   a1<-ran2[,id]
  #   max.stat <- c(max.stat,Cal_S_EqualWeight(a1)$maxS)
  # })
  
  pval <- sum(max.stat > stat1)/length(max.stat)
  return(list("pval" = pval,"max.stat" = max.stat))
}
