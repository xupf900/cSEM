setwd(" ")
library(tmvtnorm)
library(pcalg)
library(Matrix)
library(graph)
library(bnlearn)
library(BiDAG) 
library(cglasso)
library(psych)
source("sub.functions.R")
insertSource("spacefns.R",package = "BiDAG")
insertSource("usrscorefns_censoredG.R",package = "BiDAG")
insertSource("initpar.R",package = "BiDAG")
insertSource("scoreagainstdag.R",package = "BiDAG")
source("ordinalScore_censoredG.R")

###MetabolomicsData###
data<-read.csv(file='MetabolomicsData.csv')

###data sacle###
XX_sqrt <- scale(sqrt(data1))

###########################################
prob_cens=0.1      #censored percentage c=0.1,0.2,0.3
up <- quantile(XX, probs = 1 - prob_cens) # right censoring values
rr<-matrix(NA,nrow(XX),ncol(XX))
cX<-mX<-XX

for(i in 1:nrow(XX)){
  for(j in 1:ncol(XX))
    if(XX[i,j]>=up){
      rr[i,j]<-1
      cX[i,j]<- up
      mX[i,j] <- NA
    }
  else{
    rr[i,j]<-0.0
  }
}

lower=rep(-Inf,ncol(XX))
upper=rep(up,ncol(XX))

#####Calculate  S2 for complete data###########
set.seed(200)
param <- scoreparameters("usr",data=XX,usrpar = list(penType = "EBIC",
                                                     gamma = 0.5,
                                                     preLevels = NULL)); # cat("names(param) = \n");print(names(param))
param$hidden_data <- XX
param$Sigma_hat   <- cov(XX)
if (pp >= 30) {
  nr_plus1it <- 2
} else {
  nr_plus1it <- 10
}

currentDAGobj <- iterativeMCMC(param, plus1it = nr_plus1it, hardlimit = 12, alpha = 0.05)

imcmc_adj <- as.matrix(currentDAGobj$DAG)
tB.CPDAG<- dag2cpdagAdj((imcmc_adj!=0)*1)
re.mle <- PrecisionMatrix.mle(A=imcmc_adj, S=cov(XX))
S2<-S_com<- solve(re.mle$H)


###########cSEM############
##KK=20,50,100,200
source("sub.functions.R")
#set.seed(200)
cgbn<-cGBN.MCEM(cX, rr, lower, upper, max.alg = "iterativeMCMC", gamma=0.5, initial=NULL, inital.alg = "zero.matrix", KK=20, thr = 1e-02, iter.max=50,iterMCMC_alpha = 0.05, 
                usrpar = list(penType = "EBIC",
                              gamma = 0.5,
                              preLevels = NULL))
eB.CPDAG<- dag2cpdagAdj((cgbn$B!=0)*1)
(shd<-shd_cpdag(tB.CPDAG, eB.CPDAG))

S1<-S_cgbn<-solve(cgbn$H)
(S_mse<-mean((S1 -S2 )^2))

########missforest####
library("missForest")
sum<-0
k=20
for (i in 1:k) {
  out_missForest <- missForest(mX)
  X_missforest<-out_missForest$ximp
  sum=sum+X_missforest
}
X.com<-sum/k
param <- scoreparameters("usr",data=X.com,usrpar = list(penType = "EBIC",
                                                        gamma = 0.5,
                                                        preLevels = NULL)); # cat("names(param) = \n");print(names(param))
param$hidden_data <- X.com
param$Sigma_hat   <- cov(X.com)
if (pp >= 30) {
  nr_plus1it <- 2
} else {
  nr_plus1it <- 10
}

currentDAGobj <- iterativeMCMC(param, plus1it = nr_plus1it, hardlimit = 12, alpha = 0.05)
missforest_adj <- as.matrix(currentDAGobj$DAG)
eB.CPDAG<- dag2cpdagAdj((missforest_adj!=0)*1)
(shd<-shd_cpdag(tB.CPDAG, eB.CPDAG))

re.mle <- PrecisionMatrix.mle(A=missforest_adj, S=cov(X.com))
S1<-missforest_sigma<- solve(re.mle$H)
(S_mse<-mean((S1 -S2 )^2))

############missknn#####################
library("VIM")
out_knn <- kNN(mX)
X.com<- as.matrix(out_knn[, 1:pp])
param <- scoreparameters("usr",data=X.com,usrpar = list(penType = "EBIC",
                                                        gamma = 0.5,
                                                        preLevels = NULL)); # cat("names(param) = \n");print(names(param))
param$hidden_data <- X.com
param$Sigma_hat   <- cov(X.com)
if (pp >= 30) {
  nr_plus1it <- 2
} else {
  nr_plus1it <- 10
}
currentDAGobj <- iterativeMCMC(param, plus1it = nr_plus1it, hardlimit = 12, alpha = 0.05)
knn_adj <- as.matrix(currentDAGobj$DAG)
eB.CPDAG<- dag2cpdagAdj((knn_adj!=0)*1)
(shd<-shd_cpdag(tB.CPDAG, eB.CPDAG))

re.mle <- PrecisionMatrix.mle(A=knn_adj, S=cov(X.com))
S1<-knn_sigma<- solve(re.mle$H)
(S_mse<-mean((S1 -S2 )^2))

