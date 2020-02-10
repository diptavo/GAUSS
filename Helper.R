Cal_S_EqualWeight<-function(Z){
  
  K<-length(Z)
  names(Z)<-1:K
  Z1<-sort(Z, decreasing=TRUE)
  S<-rep(0,K)
  
  for(i in 1:K){
    S[i]<-sum(Z1[1:i])/sqrt(i)
  }
  id.max<-which(S == max(S,na.rm =T))
  re<-list(maxS = S[id.max], idx = as.numeric(names(Z1[1:id.max])), S=S)
  return(re)
}

CMS <- function(Z){
  
  K<-length(Z)
  names(Z)<-1:K
  Z1<-sort(Z, decreasing=TRUE)
  S<-rep(0,K)
  
  for(i in 1:K){
    S[i]<-sum(Z1[1:i])/sqrt(i)
  }
  id.max<-which(S == max(S,na.rm =T))
  re<- S[id.max]
  return(re)
}


Resample.PV <- function(stat1,cor.z,prec = 100000){
  
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


Calc.PV <- function(stat1,cor.z,max.prec = 1e+06,init = 1000,step = 1000,oc = 2,ct=250,ct.min = 100,stp = -10,tl = 500,thr = 0.05,verbose = T,LOGF = "out.log"){
  
  if(verbose)
    cat("Running resampling \n\n",file = LOGF,append = T)
  precs <- init; gpd.used = FALSE;r <- Resample.PV(stat1,cor.z,prec = precs)
  pval <- r$pval;max.stat <- r$max.stat
  while(pval < (2/init) && precs < max.prec){
    precs <- precs*step
    r <- Resample.PV(stat1,cor.z,prec = precs)
    pval <- r$pval; max.stat <- r$max.stat
  }
  if(pval < oc/length(max.stat)){
    pval <- gpd.pv(max.stat,tst = stat1,ct=ct,ct.min = ct.min,stp = stp,tl = tl,thr = thr,LOGF = LOGF);
    gpd.used = TRUE
  }
  
  return(list("pval" = pval,"gpd.used" = gpd.used))
}



mat.sqrt <- function (A)
{
  ei <- eigen(A)
  d <- ei$values
  d <- (d + abs(d))/2
  d2 <- sqrt(d)
  ans <- ei$vectors %*% diag(d2) %*% t(ei$vectors)
  return(ans)
}

condNormal <- function(x.given, mu, sigma, given.ind, req.ind){
  # Returns conditional mean and variance of x[req.ind] 
  # Given x[given.ind] = x.given 
  # where X is multivariate Normal with 
  # mean = mu and covariance = sigma 
  # 
  B <- sigma[req.ind, req.ind]
  C <- sigma[req.ind, given.ind, drop=FALSE]
  D <- sigma[given.ind, given.ind]
  CDinv <- C %*% solve(D)
  cMu <- c(mu[req.ind] + CDinv %*% (x.given - mu[given.ind]))
  cVar <- B - CDinv %*% t(C)
  list(condMean=cMu, condVar=cVar)
}

Resample.PV.C <- function(stat1,cor.z,condMean,prec = 10000){
  
  
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


Calc.PV.C <- function(stat1,cor.z,condMean,max.prec = 1e+05,init = 1000,step = 100,LOGF = "out.log"){
  
  cat("Running resampling",file = LOGF,append = T)
  precs <- init
  pval <- Resample.PV.C(stat1,cor.z,condMean = condMean,prec = precs)
  while(pval < (2/init) && precs < max.prec){
    precs <- precs*step
    pval <- Resample.PV(stat1,cor.z,prec = precs)
  }
  return(pval)
}

