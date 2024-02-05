setwd("C:/Users/lenovo/Desktop/cSEM")
rm(list=ls())
# required packages & Source ==============
library(pcalg) 
library(BiDAG)

insertSource("spacefns.R",package = "BiDAG")
insertSource("usrscorefns_censoredG.R",package = "BiDAG")
insertSource("initpar.R",package = "BiDAG")
insertSource("scoreagainstdag.R",package = "BiDAG")
source("ordinalScore_censoredG.R")

# general setting ==============
#set.seed(200)
nn = 200 # sample size
pp = 10 # the nubmer of variables in Gaussian distribution
dd = 3; # the expected sum of the in- and out-degree
prob_cens=0.5 # the precentage of censored data
lower = rep(0, pp)
upper = rep(Inf, pp)
loop=100

# data generating ==============
source("sub.functions.R")
gpData<-generate.true.graph.Para.Data(nn, pp, dd, method.dag="er", omb=0.5, omu=1, Bb=0.5, Bu=1,lower=rep(0,pp), upper=rep(Inf,pp))

  plot(gpData$trueG)
  B = gpData$truePara$B
  cX= gpData$cX # observed variables with censored value substituted by upper or lower 
  rr= gpData$rr # censored pattern in cX
  XX= gpData$XX # not censored data

#1.Our cSEM algorithm (The sampling sizes are chosen as KK = 20, 50, 100, 200 and the initial structures are chosen as inital.alg="zero.matrix, ges, mmhc") ==============
  source("sub.functions.R") 
  cgbn<-cSEM.MCEM(cX, rr, lower, upper, max.alg = "iterativeMCMC", gamma=0.5, initial=NULL, inital.alg = "ges", KK=20, thr = 1e-02, iter.max=50,iterMCMC_alpha = 0.05, 
                    usrpar = list(penType = "EBIC",
                                  gamma = 0.5,
                                  preLevels = NULL))
  
  tB.CPDAG <- dag2cpdagAdj((B!=0)*1)# change adjacent matrix B of a DAG to a adjcent matrix of its CPDAG 
  eB.CPDAG<- dag2cpdagAdj((cgbn$B!=0)*1)
  (shd.temp<-shd_cpdag(tB.CPDAG, eB.CPDAG))

#2.1 iterativeMCMC with ebic gamma = 0.5 is run on non-censored data set ==============
  library(BiDAG) 
  param <- scoreparameters("usr",data=XX,usrpar = list(penType = "EBIC",
                                                       gamma = 0,
                                                       preLevels = NULL)); # cat("names(param) = \n");print(names(param))
  param$hidden_data <- XX
  param$Sigma_hat   <- cov(XX)
  if (pp >= 30) {
    nr_plus1it <- 2
  } else {
    nr_plus1it <- 10
  }
  currentDAGobj <- iterativeMCMC(param, plus1it = nr_plus1it, hardlimit = 12, alpha = 0.05)
  
  tB.CPDAG <- dag2cpdagAdj((B!=0)*1)# change adjacent matrix B of a DAG to a adjcent matrix of its CPDAG 
  adj <- as.matrix(currentDAGobj$DAG)
  eB12.CPDAG<- dag2cpdagAdj((adj!=0)*1)
  (shd_iterativeMCMC_com_temp<-shd_cpdag(tB.CPDAG, eB12.CPDAG))
 
#2.2 iterativeMCMC with ebic gamma = 0.5 is run on the data set where the censored values are substituted by the threshold  ==============
  param <- scoreparameters("usr",data=cX,usrpar = list(penType = "EBIC",
                                                       gamma = 0,
                                                       preLevels = NULL)); # cat("names(param) = \n");print(names(param))
  param$hidden_data <- cX
  param$Sigma_hat   <- cov(cX)
  if (pp >= 30) {
    nr_plus1it <- 2
  } else {
    nr_plus1it <- 10
  }
  currentDAGobj <- iterativeMCMC(param, plus1it = nr_plus1it, hardlimit = 12, alpha = 0.05)
  
  tB.CPDAG <- dag2cpdagAdj((B!=0)*1)# change adjacent matrix B of a DAG to a adjcent matrix of its CPDAG 
  adj <- as.matrix(currentDAGobj$DAG)
  eB12.CPDAG<- dag2cpdagAdj((adj!=0)*1)
  (shd_iterativeMCMC_sub_temp<-shd_cpdag(tB.CPDAG, eB12.CPDAG))
  
#3.1 ges with ebic gamma = 0.5 is run on non-censored data set ==============
  gamma = 0.5
  score.ges <- new("GaussL0penObsScore", data = XX,
                   lambda = 0.5*log(nn) + 2*gamma*log(pp), #It is equal to 0.5*KK*log(nrow(X.complete)/KK),  
                   intercept = TRUE, use.cpp = TRUE)
  eG <- ges(score=score.ges, fixedGaps = NULL, # if fixedGaps[i,j]=True, then no edge between nodes i and j.
            adaptive = "vstructures", 
            phase = c("forward", "backward"), # "turning"), 
            iterate = TRUE, maxDegree = integer(0), 
            verbose = FALSE)
  adj <- as(eG$repr,"matrix")*1 # adj[i,j]=1 means i->j
  tB.CPDAG <- dag2cpdagAdj((B!=0)*1)# change adjacent matrix B of a DAG to a adjcent matrix of its CPDAG 
  eB33.CPDAG<- dag2cpdagAdj((adj!=0)*1)
  (shd_gescom_temp<-shd_cpdag(tB.CPDAG, eB33.CPDAG))
  

#3.2 ges with ebic gamma = 0.5 is run on the data set where the censored values are substituted by the threshold  ==============
  gamma = 0.5
  score.ges <- new("GaussL0penObsScore", data = cX,
                   lambda = 0.5*log(nn) + 2*gamma*log(pp), #It is equal to 0.5*KK*log(nrow(X.complete)/KK),  
                   intercept = TRUE, use.cpp = TRUE)
  eG2 <- ges(score=score.ges, fixedGaps = NULL, # if fixedGaps[i,j]=True, then no edge between nodes i and j.
             adaptive = "vstructures", 
             phase = c("forward", "backward"), # "turning"), 
             iterate = TRUE, maxDegree = integer(0), 
             verbose = FALSE)
  adj2 <- as(eG2$repr,"matrix")*1 # adj[i,j]=1 means i->j
  tB.CPDAG <- dag2cpdagAdj((B!=0)*1)# change adjacent matrix B of a DAG to a adjcent matrix of its CPDAG 
  eB2.CPDAG<- dag2cpdagAdj((adj2!=0)*1)
  (shd_gessub_temp<-shd_cpdag(tB.CPDAG, eB2.CPDAG))
  
#4.1. mmhc is run on non-censored data set ==============
   eM <- bnlearn::mmhc(x=data.frame(XX),maximize.args=list(score="ebic-g"))
  shd_mmhccom_temp<-shd_cpdag(tB.CPDAG, dag2cpdagAdj(bnlearn::amat(eM)))
  
# 4.2. mmhc is run on the data set where the censored values are substituted by the threshold  ==============
  
  eM2 <- bnlearn::mmhc(x=data.frame(cX),maximize.args=list(score="ebic-g"))
  shd_mmhcsub_temp<-shd_cpdag(tB.CPDAG, dag2cpdagAdj(bnlearn::amat(eM2)))
  
#5.1 pc alpha =0.001 is run on non-censored data set ==============
 
  suffStat <- list(C = cor(XX), n = nrow(XX))
  pc.fit <- pc(suffStat, indepTest = gaussCItest, p = ncol(XX), alpha = 0.001, m.max = Inf, u2pd = "relaxed",skel.method = "stable")
  shdpc_temp<-shd_cpdag(tB.CPDAG, as(pc.fit@graph, "matrix"))
  
#5.2 pc with alpha =0.001 is run on the data set where the censored values are substituted by the threshold ==============
  
  suffStat <- list(C = cor(cX), n = nrow(cX))
  pc.fit <- pc(suffStat, indepTest = gaussCItest, p = ncol(cX), alpha = 0.001, m.max = Inf, u2pd = "relaxed",skel.method = "stable")
  (shdpc_temp<-shd_cpdag(tB.CPDAG, as(pc.fit@graph, "matrix")))




