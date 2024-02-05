### This function returns the objects needed to evaluate the user defined score
usrscoreparameters <- function(initparam,
                               usrpar = list(penType = c("AIC","BIC","other"),
                                             L = 5,
                                             lambda = 2,
                                             pctesttype = "usr")){
  
  cat("============ usrscoreparameters ==================\n")
  
  return(cenGauss_ScoreParam(initparam, usrpar))
  
}

### This function evaluates the log score of a node given its parents
usrDAGcorescore <- function (j,parentnodes,n,param) {
  # cat("========= usrDAGcorescore ===========\n")
  # This function does not need to be changed, and it can be used for censored Gaussian data, as it depends on the sample covariance matrix
  return(cenGauss_CoreScore(j,parentnodes,n,param)) 
}

### This function defines the CI tests for the starting skeleton
usrdefinestartspace <- function(alpha,param,cpdag,n){
  # This function does not need to be changed, and it can be used for censored Gaussian data, as it depends on the sample covariance matrix
  # cat("========= usrdefinestartspace ===========\n")
  
  cormat<-cov2cor(param$Sigma_hat)
  if(max(abs(cormat))>1)stop("max(abs(cormat))>1")
  N <- param$N
  if(cpdag){
    pc.skel<-pcalg::pc(suffStat = list(C = cormat, n = N),
                       indepTest = gaussCItest,
                       alpha=alpha,labels=colnames(param$data),skel.method="stable",verbose = FALSE)
  } else {
    pc.skel<-pcalg::skeleton(suffStat = list(C = cormat, n = N),
                             indepTest = gaussCItest,
                             alpha=alpha,labels=colnames(param$data),method="stable",verbose = FALSE)
  }
  
  return(pc.skel)
}
