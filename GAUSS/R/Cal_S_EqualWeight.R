Cal_S_EqualWeight <-
function(Z){
  
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
