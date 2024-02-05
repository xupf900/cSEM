library(dplyr)
library(leaps)
library(pcalg)
library(glmnet)






##' cenGauss_ScoreParam(initparam,ordType): a function to be called by 'scoreparameters' in BiDAG
##' It prepares the parameters needed to calculate the score for cenGauss_ data
##' @param initparam: parameters returned from 'scoreparameters'
##' @param penType: type of penalization to be applied ("AIC" or "BIC" or "other")
##' @param L: number of truncated multivariate copies to be generated for each observation
##' @param lambda: coefficient for the BIC penalty term (default: 0 no penalty)
##' @return the input parameters
##' + estimated thresholds
##' + estimated latent Gaussian correlation matrix
cenGauss_ScoreParam <- function(initparam,
                              usrpar = list(penType = c("AIC","BIC","other"),
                                            gamma = 1,
                                            preLevels = NULL)) {

  n <- initparam$n
  initparam$penType <- usrpar$penType
  initparam$N <- nrow(initparam$data)
  initparam$gamma <- usrpar$gamma

  # Convert to matrix/array if necessary
  if (!is.matrix(initparam$data)) {
    if (is.data.frame(initparam$data)) {
      initparam$data <- as.matrix(initparam$data)
    } else {
      stop("cenGauss_ data table is not a matrix")
    }
  }
  


  initparam$hidden_data <- initparam$data
  initparam$Sigma_hat <- diag(initparam$n)


  return(initparam)
}

##' ordinalCoreScore(j,parentnodes,n,param): a function to be called by 'DAGcorescore' in BiDAG
##' It computes the ordinal score of a node given its parents
##' @param j: index of node j
##' @param parentnodes: indices of parent nodes of j
##' @param n: dimension of the variable
##' @param param: parameters returned from 'scoreparameters'
##' @return the cenGauss_ score of a node given its parents
cenGauss_CoreScore <- function(j,parentnodes,n,param) {
  # This function does not need to be changed, and it can be used for censored Gaussian data, as it depends on the sample covariance matrix

  # cat("===============cenGauss_CoreScore j =",j, "   parentnodes =", parentnodes, "    n=",  n  ,"\n")
  # print(names(param))
  if (j %in% parentnodes) {
    stop("The parent set contains the current node.")
  }

  lp <- length(parentnodes) # number of parents
  corescore <- 0
  N <- param$N

  switch(as.character(lp),
         "0" = { # no parents
           corescore <- -N/2 * (1 + log(2 * pi) + log(param$Sigma_hat[j,j])) #corescore <- -N/2 * (1 + log(2 * pi)) because in ordinal data paper, Sigma_hat[j,j]=1
         },
         "1" = { # one parent
           corescore <- -N/2 * (log(param$Sigma_hat[j,j] - param$Sigma_hat[j,parentnodes]^2 / param$Sigma_hat[parentnodes,parentnodes]) + 1 + log(2 * pi))
         },
         { # more parents
           b <- param$Sigma_hat[j,parentnodes]
           S_p <- param$Sigma_hat[parentnodes,parentnodes]
           choltemp <- chol(S_p)
           corescore <- -N/2 * (log(param$Sigma_hat[j,j] - sum(backsolve(choltemp,b,transpose=TRUE)^2)) + 1 + log(2 * pi))
         }
  )

  # cat("corescore = ", corescore, "\n")
  if (param$penType == "AIC") {
    # param$lambda <- 2 / log(param$N)
    # score <- corescore - param$lambda * log(param$N) / 2 * lp
    score <- corescore - lp
  } else if (param$penType == "EBIC") {
    # ebic = -2*loglik + log(sample.size)*|E| + 4*gamma*log(n)*|E|, where n is the number of vertices and 0 <= gamma <=1
    score <- corescore - log(param$N) / 2 * lp - 2 * param$gamma * lp * log(param$n) 
  } #else if (param$penType == "BIC") {
  #   # param$lambda <- 1
  #   # score <- corescore - param$lambda * log(param$N) / 2 * lp
  #   score <- corescore - log(param$N) / 2 * lp
  # } 
  
  # cat("param$n =", param$n, "  param$gamma=", param$gamma, "\n")
  # cat("score = ", score, "\n")

  return(score)
}


######### Taken from the source code of the pcalg package
## this function takes as parameter the adjacency matrix of a pdag (or cpdag)
## and returns the pattern of this pdag in the Meek sense, that is,
## it returns the adjacency matrix of the graph with the same skeleton where the only oriented
## edges are the v-structures (can be easily modified to work for MAGs/PAGs)
getPattern <- function(amat){

  ## makes the whole graph undirected
  tmp <- amat + t(amat)
  tmp[tmp == 2] <- 1

  ## find all v-structures i -> k <- j s.t. i not adj to k
  ## and make only those edges directed
  for (i in 1: (length(tmp[1,])-1)) {
    for (j in (i+1): length(tmp[1,])){
      if ((amat[j,i] ==0) & (amat[i,j] ==0) & (i!=j)){ ## if i no adjacent with j in G

        possible.k <- which(amat[i,]!= 0 & amat[,i]==0) ## finds all k such that i -> k is in G

        if (length(possible.k)!=0){    ## if there are any such k's then check whether j -> k for any of them
          for (k in 1: length(possible.k)){
            if ((amat[j,possible.k[k]] ==1) & (amat[possible.k[k],j]==0)) { ## if j -> k add the v-struc orientation to tmp
              tmp[possible.k[k],i] <- 0
              tmp[possible.k[k],j] <- 0
            }
          }
        }
      }
    }
  }
  tmp
}

convert_to_skeleton <- function(DAG) {
  DAG <- DAG + t(DAG)
  return(DAG != 0)
}



##' getSHD(estDAG, trueDAG):
##' a function that computes the structural difference between the PDAGs of two DAGs
##' @param estDAG: adjacency matrix of the estimated DAG
##' @param trueDAG: adjacency matrix of the true DAG
##' @return structural difference between two PDAGs
getSHD <- function(estDAG, trueDAG) {

  if (!is.matrix(estDAG)) {
    estDAG <- as(estDAG,"matrix") * 1
  }

  if (!is.matrix(trueDAG)) {
    trueDAG <- as(trueDAG,"matrix") * 1
  }

  DAGdiff <- dag2cpdag(estDAG) != dag2cpdag(trueDAG)
  return(sum(as.logical(DAGdiff + t(DAGdiff)))/2)
}


