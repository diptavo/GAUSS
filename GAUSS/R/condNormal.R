condNormal <-
function(x.given, mu, sigma, given.ind, req.ind){
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
