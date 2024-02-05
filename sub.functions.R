
dyn.load("my_shd.dll")

shd_cpdag<-function(G, H){
  #H is the PDAG by learning algorithm
  #G is the true PDAG representing equivalence
  p = ncol(G);
  if(ncol(G)!=ncol(H) ||nrow(G)!=nrow(H) ||nrow(G)!=ncol(H) )
    stop("G and H must be two square matrix with the same number of rows\n");
  dd<-.Call("Structural_Hamming_Distance", as.double(H), as.double(G), as.integer(p) );
  d<-matrix(dd, ncol=6, byrow=TRUE);
  colnames(d) <- c("SHD", "me", "ee", "md", "ed", "rd");
  return(d)
}

dag2cpdagAdj <- function(Adj){
  d <- as(Adj, "graphNEL")
  cpd <- dag2cpdag(d)
  result <- as(cpd, "matrix")
  return(result)
}

generate.true.graph.Para.Data<-function(nn, pp, dd, method.dag="er",
                                  omb=0.5,omu=1, Bb=0.5, Bu=1,
                                  lower=rep(-2,pp), 
                                  upper=rep(2,pp)
                                  ){
  
  
  # =======================================================
  # 1. generate true parameter values
  # =======================================================
  # pp  # the nubmer of variables in Gaussian distribution
  # dd  # the expected sum of the in- and out-degree
  #============== 1.1 generate DAG randomly
  trueDAG = pcalg::randDAG(n=pp, d=dd, method =method.dag, par1=NULL, par2=NULL,
                           DAG = TRUE, weighted = TRUE, wFUN = list(runif, min=0.1, max=1))
 
  B <- as.matrix(graph::graph2SparseM(trueDAG))
  
  #============== 1.2 generate parameters randomly
  mu = rep(0, pp) # mean
  om = runif(n=pp, min=omb, max=omu) # variance of noise
  beta   = runif(n=sum(B), min=Bb, max=Bu)
  pn1 = sample(x=c(1,-1), size=sum(B), replace = T)
  B[B==1] <- beta*pn1 # B[i,j] is draw uniformly from (0.5, 1) and (-1, -0.5)
  # View(B)
  Id <- diag(rep(1,pp)) # identity matrix
  Sigma <- solve(t(Id-B)) %*% diag(om) %*% solve(Id - B) #cov(X)
  tau<-qnorm(prob_cens)
  lower=(sqrt(diag(Sigma))*tau)
  HH =  (Id-B)%*%diag(1/om)%*%t(Id-B) #solve(cov(X))
  
  # =======================================================
  # 2. draw samples from censored Gaussian BN
  # =======================================================
  # =========== 2.1 draw samples from multivariate normal distribution
  Xmu = solve(Id - t(B))%*%as.matrix(mu, ncol=1) # mean of X
  XX = MASS::mvrnorm(n=nn, mu=Xmu, Sigma = Sigma)
  colnames(XX) = paste("V", 1:pp, sep="")
  
  # =========== 2.2 draw samples from censored multivariate normal distribution
  cX = XX
  rr = matrix(0, nrow=nn, ncol=pp)
  for(i in 1:pp){
    # censored below
    index = which(XX[,i]<lower[i])
    cX[index, i] = lower[i]
    rr[index, i] = -1
    
    # censored above
    index = which(XX[,i]>upper[i])
    cX[index, i] = upper[i]
    rr[index, i] = 1
  }
  
  # compute stat for cX 
  sum(rr==1)/(nn*pp)
  sum(rr==0)/(nn*pp)
  sum(rr==-1)/(nn*pp)
  
  co = rowSums(abs(rr))
  length(which(co>0))/nn
  
  truePara = list(mu=mu, B=B, om=om, Sigma=Sigma, H=HH, 
                  lower=lower, upper=upper)
  list(cX=cX, rr=rr, XX=XX, 
       trueG=trueDAG, truePara=truePara)
}

# Expectation–maximization algorithm for mle for censored univariate normal distribution within an interval.
em.unorm.estimate<-function(cX, rr, lower, upper, KK=2000, threshold=1e-2, iteration.Max = 50){
  # KK is the number of samples drawn in MCEM
  # # =========== 1. initial value ===========
  mu <- sigma <- loglik <- rep(0, iteration.Max)
  sigma<- rep(0, iteration.Max)
  # initial value in EM
  mu[1]    = mean(cX)
  sigma[1] = var(cX)
  
  # =========== 2. previous computation ====
  (n = length(cX)) # sample size
  x <- cX[rr==0] # observed samples
  nL = sum(rr==-1) # number of censored below samples
  nU = sum(rr== 1) # number of censored upper samples
 
  so = sum(x)   #sum of observed samples
  s2o= sum(x^2) #sum square of observed samples
  
  # =========== 3. EM iteration ============
  # Expectation–maximization algorithm
  
  for(iter in 1:iteration.Max){
    # E-step
    if(nL>0){
      xL <- tmvtnsim::rtnorm(mean=mu[iter], sd = sqrt(sigma[iter]), lower=-Inf, upper=lower, n = nL*KK)
      sL <- sum(xL)/KK
      s2L<- sum(xL^2)/KK
      # HL = 0 # value of H function
      Z = pnorm((lower-mu[iter])/sqrt(sigma[iter]))
      HL = nL*log(Z) + 0.5*nL*log(2*pi) + 0.5*nL*log(sigma[iter]) + mean((xL-mu[iter])^2)*nL/(2*sigma[iter])
    }
    else{
      sL <- s2L <- 0
      HL = 0 # value of H function
    }
      
    if(nU>0){
      xU <- tmvtnsim::rtnorm(mean=mu[iter], sd = sqrt(sigma[iter]), lower=upper, upper=Inf,  n = nU*KK)
      sU <- sum(xU)/KK
      s2U<- sum(xU^2)/KK
      # HU = 0 # value of H function
      Z = 1 - pnorm((upper-mu[iter])/sqrt(sigma[iter]))
      HU = nU*log(Z) + 0.5*nU*log(2*pi) + 0.5*nU*log(sigma[iter]) + mean((xU-mu[iter])^2)*nL/(2*sigma[iter])
    }
    else{
      sU <- s2U <- 0
      HU = 0 # value of H function
    }
    
    # compute log-likelihood
    fenzi = (n*mu[iter]^2) -2*mu[iter]*(so+sL+sU) + (s2o+s2L+s2U)
    
    Q <- -0.5*log(2*pi) - n/2*log(sigma[iter]) - fenzi/(2*sigma[iter])
    H <- HL + HU
    loglik[iter] <- Q + H
    
    # M-step
    
    mu[iter+1] = (so + sL + sU)/n
    sigma[iter+1] = (s2o + s2L + s2U - n*mu[iter+1]^2)/n
    
    if(abs(mu[iter+1]-mu[iter])<1e-2 && abs(sigma[iter+1]-sigma[iter])<1e-2){
      break
    }
  }
  list(mu=mu[iter+1], sigma=sigma[iter+1], 
       mu.list=mu[1:(iter+1)], sigma.list=sigma[1:(iter+1)], 
       iter=iter, loglik=loglik[1:iter])
  
}

PrecisionMatrix.mle <- function(A, S){
  # S is sample covariance matrix
  p = nrow(S)
  D<-rep(0, p)
  for(i in 1:p){
    pa<-which(A[,i]!=0)
    # cat("_____ i =",i, "   pa = ", pa, "\n")
    if(length(pa)>0){
      # fit<-lm(x[,i]~x[,pa])
  
      A[pa, i]<- solve(a= S[pa, pa], b=S[pa, i])
      D[i] <- S[i,i] - S[i, pa]%*%A[pa,i]
    }
    else{
      D[i] <- S[i,i]
    }
  }
 
  
  Id <- diag(rep(1, p))
  H <- (Id - A)%*%diag(1/D)%*%t(Id - A) # H is the precision matrix of Gaussian distribution
  
  list(H=H, A=A, D=D)
}

update.parameter.by.mle <- function(KK, X.complete, adj, S.cov, gamma){
  pp = ncol(X.complete)
  nK = nrow(X.complete)
  
  # based on sample covariance
  re.mle <- PrecisionMatrix.mle(A=adj, S=S.cov)
  H <- re.mle$H
  B <- re.mle$A
  om<- re.mle$D
  Xm<- colMeans(X.complete)
  Em<- (diag(1,ncol(X.complete)) - t(B))%*%matrix(Xm, ncol=1)# Em = (I - t(B))%*%Xm
  
  loglik <- compute.loglik(n=nK/KK, S=S.cov, wi=H)
  df = sum(B!=0) # the number of edges in DAG
  ebic <- -2*loglik + log(nn)*df + 4*gamma*df*log(pp)
  
  
  return(list(Xm=Xm, Em=Em, B=B, om=om, H=H, ebic=ebic))
}

# Compute log likelihood of Multivariate normal distribution
compute.loglik<-function(n, S, wi){
  # n is sample size
  # S is sample covariance matrix
  # wi is the inverse covariance matrix
  -0.5*n*log(2*pi) + 0.5*n*log(det(wi)) - 0.5*n*sum(diag(S%*%wi))
}

# compute conditional mean, covariance matrix, precision matrix
cond.mean.Precision <-function(mu, H, o, c, cX){
 
  if(length(o)>0){
    
    if(length(c) == 1){
      muc = mu[c] - solve(H[c,c])%*% H[c,o]%*%matrix(cX[o] - mu[o], ncol=1)
    }
    else{
      muc = mu[c] - solve(H[c,c], H[c,o])%*%matrix(cX[o] - mu[o], ncol=1)
    }
   
    Hcc = H[c,c]
    # Scc = solve(Hcc)
    
    # cat("class(Hcc)=\n");print(class(Hcc))
    # cat("Hcc=\n");print(round(Hcc,3))
    # cat("muc =", muc,"\n")
    
    if(length(c)>1){
      Hcc[abs(Hcc)<1e-7]<-0
      Hcc = as(Hcc, Class="dgCMatrix")# sparseMatrix(i = ep, j = ep, x=, dims=c(,))
    }
  }
  else{
    if(length(c) == length(mu)){
      muc = mu
      Hcc = H
      Hcc[abs(Hcc)<1e-7]<-0
      Hcc = as(Hcc, Class="dgCMatrix")# sparseMatrix(i = ep, j = ep, x=, dims=c(,))
    }
    else
      stop("length(o) = 0, but length(c) != length(mu)")
    
  }
  
  return(list(muc=muc, Hcc=Hcc))
}

LcUc<-function(lower, upper, rr, c){
  nnc = length(c)# number of censored variables
  Lc <- Uc<- rep(0, nnc)
  j <- 0
  for(i in c){
    j = j+1
    if(rr[i]==1){
      Lc[j] = upper[i]
      Uc[j] = Inf
    }
    else if(rr[i]==-1){
      Lc[j] = -Inf
      Uc[j] = lower[i]
    }
  }
  list(lower=Lc, upper=Uc)
}

mcE_step <- function(X.complete, lower, upper, nn, rr, KK, iter, seq.Xm, seq.H, cX, cputime){
  # draw samples for conditional distribution to get complete data set
  for(i in 1:nn){
    proc30 = proc.time()[3]
    o = which(rr[i,]==0)
    c = which(rr[i,]!=0)
    ind = seq(from = i, to=KK*nn, by=nn)
    # ==== 1. draw from truncated Univariate Normal Distribution
    if(length(c)==1){
      mp = cond.mean.Precision(mu=seq.Xm[iter,], H=seq.H[[iter]], o=o, c=c, cX=cX[i,])
      
      if(rr[i,c]==-1)
        xs <- tmvtnsim::rtnorm(mean=mp$muc, sd = sqrt(1/mp$Hcc), lower=-Inf,     upper=lower[c], n = KK)
      if(rr[i,c]==1)
        xs <- tmvtnsim::rtnorm(mean=mp$muc, sd = sqrt(1/mp$Hcc), lower=upper[c], upper=Inf,      n = KK)
      X.complete[ind,c] <- xs
     
    }
    proc31 = proc.time()[3]
    cputime["E1"] = cputime["E1"] + proc31 - proc30
    
    # ==== 2. draw from truncated Multivariate Normal Distribution
    if(length(c)>=2){
      LU = LcUc(lower, upper, rr[i,], c); # cat("LU = \n"); print(LU)
      mp = cond.mean.Precision(mu=seq.Xm[iter,], H=seq.H[[iter]], o=o, c=c, cX=cX[i,])
      
      time0 = proc.time()[3]
      # draw samples from truncated Multivariate Normal Distribution
      xs <- tmvtnorm::rtmvnorm.sparseMatrix(n=KK, mean = mp$muc, H=mp$Hcc, 
                                            lower=LU$lower, upper=LU$upper, burn.in.samples=1000, thinning=10)
      X.complete[ind,c] <- xs
      time1 = proc.time()[3]
      cputime["E2"] = cputime["E2"] + time1 - time0
    }
    
  }
  return(list(X.complete=X.complete, cputime=cputime))
}

cSEM.MCEM<-function(cX, rr, lower, upper, max.alg = iterativeMCMC, gamma =0.5, initial=NULL, 
                    inital.alg = c("zero.matrix", "mmhc", "ges"), KK=500, thr = 1e-03, iter.max=50,
                    iterMCMC_alpha = 0.05, hardlimit = 12,
                    usrpar = list(penType = c("AIC", "EBIC"),
                                  gamma = 0.5,
                                  preLevels = NULL)){
  pp = ncol(cX)
  nn = nrow(cX)
  
  if(max.alg !="iterativeMCMC")
    stop("max.alg !=iterativeMCMC")

  if(max.alg == "iterativeMCMC"){
    param <- scoreparameters("usr",data=cX,usrpar = usrpar); # cat("names(param) = \n");print(names(param))
    
    # Maximum iterations to control runtime (can change)
    if (pp >= 30) {
      nr_plus1it <- 2
    } else {
      nr_plus1it <- 10
    }
    gamma = usrpar$gamma
  }
  if(gamma>1 || gamma<0)stop("In ebic, 0 <= gamma <= 1.")
  
  # ====================== CPU time =====================
  cputime <-rep(0, 12)
  names(cputime)<-c("initial", "E1", "E2", "E3", "E-step", "M.ges.score", "M.ges", "M.mmhc", "M.iterMCMC", "M.mle", "M-step", "total")
  
  proc1 = proc.time()[3]
  # ====================== intermediate parameter =======
  seq.B <- seq.H <- seq.S <- list()
  seq.Xm <- seq.Em <- seq.om <- matrix(0, nrow=iter.max, ncol=pp)
  ebic <- rep(Inf, iter.max)
 
  # X.complete is used to store the complete samples in E-step
  X.complete = cX 
  if(KK >= 2){
    for(i in 2:KK){
      X.complete = rbind(X.complete, cX)
    }
  }
  X2.complete = X.complete
  nK = nrow(X.complete)
  
  # ====================== compute initial values ===============
  # == compute an initial B 
  if(is.null(initial)){
    if(inital.alg == "zero.matrix"){
      B = matrix(0, nrow=pp, ncol=pp)
    }
    else if (inital.alg == "mmhc"){
      ## ==== compute an initial B by mmhc according to BIC score =====
      eG <- bnlearn::mmhc(x=data.frame(cX)) 
      adj <- bnlearn::amat(eG) # adj[i,j]=1 means i->j
      re.mle <- PrecisionMatrix.mle(A=adj, S=cov(cX))
      B <- re.mle$A
    }
    else if (inital.alg == "ges"){
      ## ==== compute an initial B by ges according to BIC score =====
      score.ges <- new("GaussL0penObsScore", data = cX,
                       lambda = 0.5*log(nrow(cX)), 
                       intercept = TRUE, use.cpp = TRUE)
      
      eG <- ges(score=score.ges, fixedGaps = NULL, # if fixedGaps[i,j]=True, then no edge between nodes i and j.
                adaptive = "vstructures", 
                phase = c("forward", "backward"), # "turning"), 
                iterate = TRUE, maxDegree = integer(0), 
                verbose = FALSE)
      adj <- as(eG$repr,"matrix")*1 # adj[i,j]=1 means i->j
      re.mle <- PrecisionMatrix.mle(A=adj, S=cov(cX))
      B <- re.mle$A
    }
    else{
      stop("inital.alg = 'zero.matrix', 'mmhc', 'ges'.")
    }
    
    # compute an intial mu and omega 
    mu <- om <- rep(NA, pp)
    if(KK==1)
      KK.em=50
    else
      KK.em = KK
    for(i in 1:pp){
      result= em.unorm.estimate(cX=cX[,i], rr=rr[,i], KK=KK.em, lower=lower[i], upper=upper[i])
      mu[i] = result$mu
      om[i] = result$sigma
    }
    
    initial<-list(mu=mu, B=B, om=om, alg=inital.alg)
  }
  
  # ==== End computing initial values ===============
  
  Id <- diag(1, pp)
  seq.Xm[1,]  <- initial$mu
  seq.om[1,]  <- initial$om
  seq.B[[1]]  <- initial$B
  seq.H[[1]]  <- (Id-seq.B[[1]])%*%diag(1/seq.om[1,])%*%t(Id-seq.B[[1]])
  seq.S[[1]]  <- solve(t(Id-seq.B[[1]]))%*%diag(seq.om[1,])%*%solve(Id-seq.B[[1]])
  # record the error of two X.complete and X2.complete
  
  
  df = sum(seq.B[[1]]!=0) # the number of edges in DAG
  ebic[1] <- -2*compute.loglik(n=nn, S=cov(cX)*(nn-1)/nn, wi=seq.H[[1]]) + log(nn)*df + 4*gamma*df*log(pp) 
  
  
  # ====================== em iteration ===============
  proc2 = proc.time()[3]
  
  cputime["initial"] = proc2 - proc1
  
  for(iter in 1:iter.max){
    proc3 = proc.time()[3]
    #========== E-step: we use MC method to approximate conditional expectation
    
    result <- mcE_step(X.complete, lower, upper, nn, rr, KK, iter, seq.Xm, seq.H, cX, cputime)
    cputime    <- result$cputime
    X.complete <-result$X.complete
    S.cov <- ((nK-1)/nK)*cov(X.complete)#S.cov = sum_i=1^n sum_k=1^K(xik - xb)(xik - xb)^T/(n*K)
    
   
    proc4 = proc.time()[3]
    cputime["E-step"] = cputime["E-step"] + proc4 - proc3
    
    #========== M-step
    
    # call iterativeMCMC to maximize EBIC value in M-step
    
    if(max.alg == "iterativeMCMC"){
      proc5 = proc.time()[3]
      param$hidden_data <- X.complete
      param$Sigma_hat   <- S.cov
      
      currentDAGobj <- iterativeMCMC(param, plus1it = nr_plus1it, hardlimit = hardlimit, alpha = iterMCMC_alpha)
      adj <- as.matrix(currentDAGobj$DAG)
      
      proc6 = proc.time()[3]
      cputime["M.iterMCMC"] = cputime["M.iterMCMC"] + proc6 - proc5
    }
    
    # update parameter by maxmimum likelhood estimation and compute a EBIC value
    re <- update.parameter.by.mle(KK, X.complete, adj, S.cov, gamma)
    seq.Xm[iter+1,] = re$Xm 
    seq.Em[iter+1,] = re$Em
    seq.om[iter+1,] = re$om 
    seq.B[[iter+1]] = re$B
    seq.H[[iter+1]] = re$H
    ebic[iter+1]    = re$ebic
    
    proc7 = proc.time()[3]
    cputime["M.mle"] = cputime["M.mle"] + proc7 - proc6
    cputime["M-step"]     = cputime["M-step"]     + proc7 - proc4
 
    if(abs((ebic[iter+1]-ebic[iter])/ebic[iter]) < thr){
      break
    }
  }
  
  # ==================
  seq = list(Xm=seq.Xm[1:(iter+1),], 
             Em=seq.Em[1:(iter+1),], 
             om=seq.om[1:(iter+1),], 
             B=seq.B, H=seq.H, 
             ebic=ebic[1:(iter+1)])
  proc8 = proc.time()[3]
  cputime["total"] = proc8 - proc1
  return(list(Xmu=re$Xm, EMu=re$Em, om=re$om, B=re$B, H=re$H, 
              max.alg = max.alg,
              initial=initial, ebic=ebic[iter+1], gamma.ebic = gamma, iter=iter, seq=seq, cputime=cputime,X.complete=X.complete))
}

