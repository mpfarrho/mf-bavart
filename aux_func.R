
# Horseshoe prior in its auxiliary representation, see Makalic and Schmidt (2015)
get.hs <- function(bdraw,lambda.hs,nu.hs,tau.hs,zeta.hs){
  k <- length(bdraw)
  if (is.na(tau.hs)){
    tau.hs <- 1   
  }else{
    tau.hs <- invgamma::rinvgamma(1,shape=(k+1)/2,rate=1/zeta.hs+sum(bdraw^2/lambda.hs)/2) 
  }
  
  lambda.hs <- invgamma::rinvgamma(k,shape=1,rate=1/nu.hs+bdraw^2/(2*tau.hs))
  
  nu.hs <- invgamma::rinvgamma(k,shape=1,rate=1+1/lambda.hs)
  zeta.hs <- invgamma::rinvgamma(1,shape=1,rate=1+1/tau.hs)
  
  ret <- list("psi"=(lambda.hs*tau.hs),"lambda"=lambda.hs,"tau"=tau.hs,"nu"=nu.hs,"zeta"=zeta.hs)
  return(ret)
}

# BART regression
BART.reg <- function(Y, X, nr = 1, nsave, nburn, cgm.exp = 2, cgm.level = 0.95, num.trees=200){
  ntot <- nsave+nburn
  
  #Prepare additional matrices
  y <- Y[, nr]
  #Creates a Matrix
  if (nr > 1) Z <- Y[, 1:(nr-1), drop =F] else Z <- NULL
  
  
  N <- nrow(X)
  M <- ncol(Y)
  priormu <- c(0, 10)
  priorphi <- c(20, 1.5)
  priorsigma <- 1
  svdraw <- list(para = c(mu = -10, phi = 0.9, sigma = 0.2), latent = rep(-10, N))
  
  #Initialize HS prior on covariances (if any)
  if (nr > 1){
    lambda.A <- 1
    nu.A <- 1
    tau.A <- 1
    zeta.A <- 1
    prior.cov <- rep(1, ncol(Z))
  }
  control <- dbartsControl(verbose = FALSE, keepTrainingFits = TRUE, useQuantiles = FALSE,
                           keepTrees = TRUE, n.samples = ntot,
                           n.cuts = 100L, n.burn = nburn, n.trees = num.trees, n.chains = 1,
                           n.threads = 1, n.thin = 1L, printEvery = 1,
                           printCutoffs = 0L, rngKind = "default", rngNormalKind = "default",
                           updateState = FALSE)
  
  sampler <- dbarts(y~X, control = control,tree.prior = cgm(cgm.exp, cgm.level), node.prior = normal(10),n.samples = nsave, weights=rep(1,N))
  
  varcount <- matrix(NA_integer_, nsave, ncol(X))
  #Store some objects
  h.store <- array(NA, c(nsave, N))
  sigma2.store <- array(NA, c(nsave, 1))
  if (nr > 1) A0.store <- array(NA, c(nsave, nr - 1)) else A0.store <-  NA
  
  fit.store <- array(NA, c(nsave, N))
  start <- Sys.time()
  
  for (irep in seq_len(ntot)){
    #Step I: Simulate the constant error variance and the BART part
    rep.i <- sampler$run(0L, 1L) #Samples trees and sigma
    #plot(y, main=irep); lines(rep.i$train, col="red")
    if (any(is.na(rep.i$train))) break
    #Step II: Simulate the stochastic volatility part conditional on the fit from the tree part of the model
    svdraw <- svsample2((y - rep.i$train)+1e-8, startpara = para(svdraw), startlatent = latent(svdraw), priormu = priormu,  priorphi = priorphi, priorsigma = priorsigma)
    normalizer <- as.numeric(exp(-.5*latent(svdraw)))
    #sampler$setWeights(exp(-latent(svdraw)))
    weights.new <- as.numeric(exp(-latent(svdraw)))
    dat <- dbartsData(formula = y~X,weights=weights.new)
    sampler$setData(dat)
    #Step III: Sample covariance parameters using a standard regression
    if (nr > 1){
      y.hat  <- (y-rep.i$train)*normalizer
      Z.hat <- Z * normalizer
      if (ncol(Z) > 1) V.prior.cov.inv <- diag(1/prior.cov) else V.prior.cov.inv <- 1/prior.cov
      V.cov <- solve(crossprod(Z.hat) + V.prior.cov.inv)
      A.cov <- V.cov %*% crossprod(Z.hat, y.hat)
      A.cov.draw <- A.cov + t(chol(V.cov)) %*% rnorm(ncol(Z))
      
      sampler$setResponse(y - Z%*%A.cov.draw)
      
      #Step III.b: Sample prior shrinkage parameters using a Horseshoe
      hs_draw <- get.hs(bdraw=A.cov.draw,lambda.hs=lambda.A,nu.hs=nu.A,tau.hs=tau.A,zeta.hs=zeta.A)
      
      lambda.A <- hs_draw$lambda
      nu.A <- hs_draw$nu
      tau.A <- hs_draw$tau
      zeta.A <- hs_draw$zeta
      prior.cov <- hs_draw$psi
    }
    if (irep == nburn){
    #  control@updateState <- TRUE
    #  sampler <- dbarts(y~X, control = control,tree.prior = cgm(cgm.exp, cgm.level), node.prior = normal(2), n.samples = nsave)
    }
    #Final step: Store draws of interest
    if (irep > nburn){
      varcount[irep-nburn,] <- t(rep.i$varcount)
      sigma2.store[irep-nburn,] <- rep.i$sigma^2
      h.store[irep-nburn,] <- latent(svdraw)#+log(rep.i$sigma)
      fit.store[irep-nburn,] <- rep.i$train+rnorm(N, 0, exp(.5*latent(svdraw)))
      if (nr>1) A0.store[irep-nburn,] <- A.cov.draw
    }
    
    # if (irep %% 50==0){
    #  end <- Sys.time()
    #  print(paste0("Average time necessary to produce a single draw over the last 50 draws ", round(as.numeric(end-start)/50, digits=4), " seconds"))
    # print(irep)
    # start <- Sys.time()
    # plot(y, main=irep); lines(rep.i$train, col="red")
    # }
  }
  
  res.list <- list(BART = sampler, h = h.store, A0 = A0.store, sigma2 = sigma2.store, varcount =varcount)
  return(res.list)
}

BARTVAR <- function(Y=Y, nsave, nburn, p, fhorz, imphor, cgm.exp = 2, cgm.level = 0.95, num.trees = 75,predictions=TRUE, IRFs = FALSE, sl.shock =3){
  X <- embed(Y , p+1)
  Y <- Y[(p+1):nrow(Y),, drop=F]
  X <- X[, (ncol(Y)+1):ncol(X)] #Input X
  
  #Compute a proper labeling for the X
  # X.names <- NULL
  # for (jj in 1:p) X.names <- c(X.names, paste0(colnames(Y), "(t-",jj,")"))
  # colnames(X) <- X.names
  
  
  M <- ncol(Y)
  N <-  nrow(X)
  if (IRFs){
    shock.vec <- matrix(0, M, 1)
    shock.vec[sl.shock] <- 1
  }
  
  bartvar.list <- list()
  for (i in seq_len(M)){
    bartvar.list[[i]] <-  BART.reg(Y=Y, X=X, nr=i, nsave=nsave, nburn=nburn, cgm.exp = cgm.exp, cgm.level = cgm.level, num.trees=num.trees)
  }
  
  #Now recover the VAR structure of the model and perform forecasting
  A0.full <- matrix(0, M, M)
  H.full <- matrix(0, N, M)
  diag(A0.full) <- rep(1,M)
  
  Sigma.store <- array(NA, c(nsave, M, M))
  H.store <- array(NA, c(nsave, N, M))
  fit.store <- array(NA, c(nsave, N, M))
  if (IRFs) IRF.store <- array(NA, c(nsave, M, imphor)) else IRF.store <- NULL
  if (predictions) pred.store <- array(NA, c(nsave, fhorz, M)) else pred.store <- NULL
  for (i in seq_len(M)) fit.store[,,i] <- t(bartvar.list[[i]]$BART$predict(X))
  for (irep in seq_len(nsave)){
    print(irep)
    #Step 1: Recover the lower Cholesky factor of the variance-covariance matrix
    for (i in seq_len(M)){
      if (i > 1){
        A0.full[i, 1:(i-1)] <- -bartvar.list[[i]]$A0[irep,]
      }
      H.full[,i] <- exp(bartvar.list[[i]]$h[irep,]) * exp(bartvar.list[[i]]$sigma2[irep,])
      H.store[irep,,i] <- exp(bartvar.list[[i]]$h[irep,]) * exp(bartvar.list[[i]]$sigma2[irep,])
    }
    #Step 2: Compute fitted values for the full system
    A0.inv <- solve(A0.full)
    fit.store[irep,,] <- fit.store[irep,,]%*%t(A0.inv)
    Sigma.T <- A0.inv %*% diag(apply(H.full,2,mean)) %*% t(A0.inv)
    
    Sigma.store[irep,,] <-  Sigma.T
    
    Y.hat <- matrix(0, M, fhorz)
    #Do predictions
    if (predictions){
      if (p == 1){
        X.hat <- c(Y[N, ])
      }else{
        X.hat <- c(Y[N, ], X[N, 1:(M*(p-1))])  
      }
      
      tree.pred <- matrix(0, M, 1)
      for (nn in seq_len(fhorz)){
        for (j in seq_len(M)) tree.pred[j] <-   bartvar.list[[j]]$BART$predict(X.hat)[, irep]
        Y.tp1 <- A0.inv %*% tree.pred + t(chol(Sigma.T))%*%rnorm(M)
        
        if (p == 1){
          X.hat <- as.numeric(Y.tp1)
        }else{
          X.hat <- c(Y.tp1, X.hat[1:(M*(p-1))])
        }
        
        Y.hat[,nn] <- Y.tp1
      }
      pred.store[irep, , ] <-  t(Y.hat)
    }
    if (IRFs){
      
      if (p == 1){
        X.hat.0 <- c(Y[N, ])
        X.hat.1 <- c(Y[N, ])
      }else{
        X.hat.0 <- c(Y[N, ], X[N, 1:(M*(p-1))])
        X.hat.1 <- c(Y[N, ], X[N, 1:(M*(p-1))])
      }
      IRF <- array(NA, c(M, imphor))
      Y.new <- Y[N,]
      tree.pred.0 <- tree.pred.1 <- matrix(0, M, 1)
      for (nn in seq_len(imphor)){
        for (j in seq_len(M)){
          tree.pred.1[j] <-   bartvar.list[[j]]$BART$predict(X.hat.1)[, irep]
          tree.pred.0[j] <-   bartvar.list[[j]]$BART$predict(X.hat.0)[, irep]
        } 
        if (nn==1) shock.it <-  t(chol(Sigma.T))%*%shock.vec else shock.it <- 0
        Y.tp1 <- A0.inv %*% tree.pred.1 + shock.it
        Y.tp0 <- A0.inv %*% tree.pred.0 #All shocks zeroed out
        
        if (p == 1){
          X.hat.1 <- as.numeric(Y.tp1)
          X.hat.0 <- as.numeric(Y.tp0)
        }else{
          X.hat.1 <- c(Y.tp1, X.hat.1[1:(M*(p-1))])
          X.hat.0 <- c(Y.tp0, X.hat.1[1:(M*(p-1))])
        }
        IRF[, nn] <- Y.tp1 - Y.tp0
      }
      IRF.store[irep, ,] <- IRF
    }
  }
  
  fit <- apply(fit.store, c(2,3), mean, na.rm=TRUE)
  ts.plot(cbind(fit[,2], Y[,2]), col=c(1,2))
  
  if (IRFs){
    IRF.mean <- apply(IRF.store, c(2,3), median, na.rm=TRUE)
    IRF.low <- apply(IRF.store, c(2,3), quantile, 0.16, na.rm=TRUE)
    IRF.high <- apply(IRF.store, c(2,3), quantile, 0.84, na.rm=TRUE)
    
    H.normalized <- apply(H.store, c(1,3), function(x) (x-mean(x))/sd(x))
    H.normalized <- aperm(H.normalized, c(2,1,3))
    vola.mean <- apply(H.normalized, c(2,3), median)
    vola.low <- apply(H.normalized, c(2,3), quantile, 0.16)
    vola.high <- apply(H.normalized, c(2,3), quantile, 0.84)
    
    dir.create("Results_insample")
    nhor.max <- 20
    line.size <- 2
    line.color.bounds <- "lightgray"
    line.color.med <- "black"
    par(mfrow=c(1,1))
    for (jj in 1:M){
      #Plot IRFs
      pdf(paste0("Results_insample/IRF_", colnames(Y)[[jj]], ".pdf"))
      par(mar=c(2,2,2,2))
      ts.plot(xlab="",main="", cbind(IRF.high[jj,1:nhor.max], IRF.mean[jj, 1:nhor.max], IRF.low[jj,1:nhor.max], 0), col=c(line.color.bounds, line.color.med, line.color.bounds), lty=c(2,1,2), lwd=line.size); 
      abline(h=0, col="red", lwd=line.size)
      dev.off()
      
      #Plot SV part of the model
      pdf(paste0("Results_insample/vola_", colnames(Y)[[jj]], ".pdf"))
      par(mar=c(2,2,2,2))
      vola.mat <- ts(cbind(vola.high[,jj], vola.mean[, jj], vola.low[,jj]), start=c(1999, 7), frequency = 12)
      ts.plot(xlab="",main="", vola.mat, col=c(line.color.bounds, line.color.med, line.color.bounds), lty=c(2,1,2), lwd=line.size); 
      abline(h=0, col="red", lwd=line.size)
      dev.off()
    }
  }
  
  res.obj <- list(predictions = pred.store, IRFs = IRF.store, vola = H.store, SIGMA = Sigma.store, bart.obj = bartvar.list, fit = fit.store)
}

mix.approx <- function(Ystar, h){

#Ystar <- rnorm(100)
#h <- rnorm(100)

N <- length(Ystar)
pi <- c(0.0073, .10556, .00002, .04395, .34001, .24566, .2575)
mi <- c(-10.12999, -3.97281, -8.56686, 2.77786, .61942, 1.79518, -1.08819) - 1.2704 # %% means already adjusted!! %%
sig2 <- c(5.79596, 2.61369, 5.17950, .16735, .64009, .34023, 1.26261)


temprand <- runif(N)
z <- Ystar - h
rep(pi, N)

mi.t <- matrix(NA, N,1)
sig2.t <-matrix(NA, N,1)
for (n in seq_len(N)){
  
  prob.n <- dnorm(z[n],mi, sqrt(sig2))*pi
  prob.n <- prob.n/sum(prob.n)
  
  s.i <- sample(1:7,1, prob = prob.n)
  
  mi.t[n] <- mi[s.i]
  sig2.t[n] <- sig2[s.i]
}

return(list(mu=mi.t, sig2=sig2.t))
}

# construct the companion matrix                           
get_companion <- function(Beta_,varndxv){
  nn <- varndxv[[1]]
  nd <- varndxv[[2]]
  nl <- varndxv[[3]]
  
  nkk <- nn*nl+nd
  
  Jm <- matrix(0,nkk,nn)
  Jm[1:nn,1:nn] <- diag(nn)
  
  if(nd==1){
    MM <- rbind(t(Beta_),cbind(diag((nl-1)*nn), matrix(0,(nl-1)*nn,nn+1)),c(matrix(0,1,(nn*nl)),1))
  }else{
    MM <- rbind(t(Beta_),cbind(diag((nl-1)*nn), matrix(0,(nl-1)*nn,nn)))
  }
  
  return(list(MM=MM,Jm=Jm))
}
                          
# FFBS with lags
ffbs <- function(y,MM,MMcons,Ft,Rt,Sig_big,beta11,p11,kk,t,p,m, beta0, P00){
# 

  beta_tt <- matrix(0,t,kk)
  ptt <- array(0,dim=c(t,kk,kk))
  
  beta11 <- beta0
  p11 <- P00
  
  # filtering
  for(i in 1:t){
    beta10 <- MMcons+MM%*%beta11 # forecast
    p10 <- MM%*%p11%*%t(MM)+Sig_big[,,i] # forecast error
    
    yhat <- Ft[,,i]%*%beta10 # map to observation eq.
    eta <- y[i,]-yhat # forecast error
    fe <- (Ft[,,i]%*%p10%*%t(Ft[,,i]))+Rt[,,i] # forecast error variance
    fe_inv <- solve(fe)
    
    K <- (p10%*%t(Ft[,,i]))%*%fe_inv
    beta11 <- beta10+K%*%eta
    p11 <- p10-K%*%(Ft[,,i]%*%p10)
    
   # p11[p11>10] <- 10
    ptt[i,,] <- p11
    beta_tt[i,] <- beta11
   # if (i == 344) break()
  }
  
  # smoothing (backward recursions)
  beta2 <- matrix(0,t,kk) # collect state variable draws
  bm2 <- beta2
  jv <- 1:m # companion selector
  jv1 <- seq(1,kk,by=m) # selector of lagged states
  
  p00 <- ptt[t,jv1,jv1]
  beta2[t,] <- beta_tt[t,]
  beta2[t,jv1] <- beta_tt[t,jv1]+t(chol(p00))%*%rnorm(p)
  
  q <- Sig_big[jv,jv,t]
  mu <- MMcons[jv]
  f <- MM[jv,]
  
  for (i in (t-1):1) {
    pt <- ptt[i,,]
    fpfq_inv <- solve(f%*%tcrossprod(pt,f)+q)
    bm <- beta_tt[i,]+t(pt%*%t(f)%*%fpfq_inv%*%t(beta2[i+1,jv]-mu-beta_tt[i,]%*%t(f)))
    pm <- pt-pt%*%t(f)%*%fpfq_inv%*%f%*%pt
    beta2[i,] <- bm
    beta2[i,jv1] <- bm[jv1]+t(chol(pm[jv1,jv1]))%*%rnorm(p)
    bm2[i,] <- bm
  }
 # browser()
  out <- c(rep(0,p),beta2[,1])
  return(out)
}

KF.slow <- function(y,MM,Ft,Rt,Sig_big,kk,t, beta0, P00){
  # Y00[is.na(Y00)] <- 0
  # y=(Y00);MM = MM; Ft = Ft; Rt=Rt; Sig_big = Sig_big; kk = K; t = N; beta0 = beta0; P00 = P00
  
  beta_tt <- matrix(0,t,kk)
  ptt <- array(0,dim=c(t,kk,kk))
  
  beta11 <- beta0
  p11 <- P00
  # filtering
  for(i in 1:t){
    beta10 <- MM%*%beta11 # forecast
    p10 <- MM%*%p11%*%t(MM)+Sig_big[,,i] # forecast error
    
    yhat <- Ft[,,i]%*%beta10 # map to observation eq.
    eta <- y[i,]-yhat # forecast error
    fe <- (Ft[,,i]%*%p10%*%t(Ft[,,i]))+Rt[,,i] # forecast error variance
    fe_inv <- try(solve(fe), silent=TRUE)
    if (is(fe_inv, "try-error")) fe_inv <- MASS::ginv(fe)
    
    K.gain <- (p10%*%t(Ft[,,i]))%*%fe_inv
    beta11 <- beta10+K.gain%*%eta
    p11 <- p10-K.gain%*%(Ft[,,i]%*%p10)

    ptt[i,,] <- p11
    beta_tt[i,] <- t(beta11)
  }

  #print(c(mean(beta_tt[6,seq(1, K, M)]), y[6,1]))
  
  
  out <- list(att = t(beta_tt), Ptt = aperm(ptt, c(2,3,1)))
  return(out)
}


RATE = function(X = X, f.draws = f.draws,prop.var = 1, low.rank = FALSE, rank.r = min(nrow(X),ncol(X)), nullify = NULL,snp.nms = snp.nms, cores = 1){
  
  ### Install the necessary libraries ###
  usePackage("doParallel")
  usePackage("MASS")
  usePackage("Matrix")
  usePackage("svd")
  
  ### Determine the number of Cores for Parallelization ###
  if(cores > 1){
    if(cores>detectCores()){warning("The number of cores you're setting is larger than detected cores!");cores = detectCores()}
  }
  
  ### Register those Cores ###
  registerDoParallel(cores=cores)
  
  ### First Run the Matrix Factorizations ###  
  svd_X = svd(X,rank.r); 
  dx = svd_X$d > 1e-10
  px = cumsum(svd_X$d^2/sum(svd_X$d^2)) < prop.var
  r_X = dx&px 
  u = with(svd_X,(1/d[r_X]*t(u[,r_X])))
  v = svd_X$v[,r_X]  
  
  if(low.rank==TRUE){
    # Now, calculate Sigma_star
    SigmaFhat = cov(f.draws)
    Sigma_star = u %*% SigmaFhat %*% t(u)
    
    # Now, calculate U st Lambda = U %*% t(U)
    svd_Sigma_star = svd(Sigma_star,rank.r)
    r = svd_Sigma_star$d > 1e-10
    U = t(ginv(v)) %*% with(svd_Sigma_star, t(1/sqrt(d[r])*t(u[,r])))
    
    mu = v%*%u%*%colMeans(f.draws)
  }else{
    beta.draws = t(ginv(X)%*%t(f.draws))
    V = cov(beta.draws); #V = as.matrix(nearPD(V)$mat)
    D = ginv(V)
    svd_D = svd(D)
    r = sum(svd_D$d>1e-10)
    U = with(svd_D,t(sqrt(d[1:r])*t(u[,1:r])))
    
    mu = colMeans(beta.draws)
  }
  
  ### Create Lambda ###
  Lambda = tcrossprod(U)
  
  ### Compute the Kullback-Leibler divergence (KLD) for Each Predictor ###
  int = 1:length(mu); l = nullify;
  
  if(length(l)>0){int = int[-l]}
  
  if(nrow(X) < ncol(X)){
    KLD = foreach(j = int, .combine='c')%dopar%{
      q = unique(c(j,l))
      m = abs(mu[q])
      U_Lambda_sub = qr.solve(U[-q,],Lambda[-q,q,drop=FALSE])
      kld = crossprod(U_Lambda_sub%*%m)/2
      names(kld) = snp.nms[j]
      kld
    }
  }else{
    KLD = foreach(j = int, .combine='c')%dopar%{
      q = unique(c(j,l))
      m = mu[q]
      alpha = t(Lambda[-q,q])%*%ginv(as.matrix(nearPD(Lambda[-q,-q])$mat))%*%Lambda[-q,q]
      kld = (t(m)%*%alpha%*%m)/2
      names(kld) = snp.nms[j]
      kld
    }
  }
  
  ### Compute the corresponding “RelATive cEntrality” (RATE) measure ###
  RATE = KLD/sum(KLD)
  
  ### Find the entropic deviation from a uniform distribution ###
  Delta = sum(RATE*log((length(mu)-length(nullify))*RATE))
  
  ### Calibrate Delta via the effective sample size (ESS) measures from importance sampling ###
  #(Gruber and West, 2016, 2017)
  ESS = 1/(1+Delta)*100
  
  ### Return a list of the values and results ###
  return(list("KLD"=KLD,"RATE"=RATE,"Delta"=Delta,"ESS"=ESS))
}

### Define the Package ###
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
