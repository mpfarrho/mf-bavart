mfbart <- function(Y.raw,VAR.mean="bart",p=3,nburn=5000,nsave=5000,prior.cov=0.01,cgm.level=0.95,cgm.exp=2,sd.mu=2,num.trees=250, prior.sig = c(5, 0.99), n.quarter=1, quiet=FALSE, iter.update=50, var.thrsh=3, max.count.var = 10, exact=FALSE){
  require(dbarts)
  require(stochvol)
  require(FKF)
  require(forecast)
  source("aux_func.R")
  
  ntot <- nsave+nburn
  n.monthly <- ncol(Y.raw)-n.quarter # number of monthly quantities
  
  # interpolate missings in monthly variables
  Y.months <- Y.raw[,-1]
  id.num <- apply(is.na(Y.months),2,sum)
  id.na <- which(id.num>0)
  for(j in id.na){
    yy <- Y.months[,j]
    yy[is.na(yy)] <- as.numeric(forecast(auto.arima(yy[!is.na(yy)]),h=id.num[j])$mean)
    Y.months[,j] <- yy
  }
  Y.raw[,2:ncol(Y.raw)] <- Y.months
  
  # get dimensions
  N <- nrow(Y.raw)-p
  M <- ncol(Y.raw)
  
  # Stuff necessary for the mixed frequency part of our model
  mid <- is.na(Y.raw)*1 # id for missing observations
  
  # create observation eq. quantities
  Ft <- array(0,dim=c(M,M*p,N)) # predetermined loadings
  Rt <- array(0,dim=c(M,M,N)) # measurement error variance
  
  F0 <- matrix(0, M, M*(p))
  FX <- matrix(0, M, M*p)
  
  F0[(n.quarter+1):nrow(F0),(n.quarter+1):nrow(F0)] <- diag(n.monthly)
  FX[(n.quarter+1):nrow(F0),(n.quarter+1):nrow(F0)] <- diag(n.monthly)
  for(i in 1:n.quarter){
    sl <- seq(1, p*M, M)
    if (p == 5) loadings <-  1/9*c(1, 2, 3, 2, 1)  else loadings <- c(1/3, 1/3, 1/3)
    zero.vec <- rep(0, M*p)
    zero.vec[sl] <- loadings
    F0[i,] <- zero.vec
  }
  
  mid.short <- mid[(p+1):nrow(mid),]
  for(tt in seq_len(N)){
    for (jj in seq_len(M)){
      if(mid.short[tt,jj]==1){
        Rt[jj,jj,tt] <- 10^10
      }else{
        Rt[jj,jj,tt] <- 0
      }  
    }
    
    if(mid.short[tt,1]==1){
      Ft[,,tt] <- FX
    }else{
      Ft[,,tt] <- F0
    }
  }
  
  # Interpolate first set of missings using BART (using only monthly vars)
  X.bart <- as.matrix(Y.raw[, 2, drop=F])
  class(X.bart) <- "numeric"
  Y.bart <- Y.raw[,1]
  
  for (j in seq_len(1)){
    interp.bart <- bart(X.bart,as.numeric(Y.bart),keeptrees=TRUE)
    fit.j <- apply(predict(interp.bart, newdata=X.bart), 2, mean)
    
    Y.raw[mid[, j]==1, j] <- fit.j[mid[ , j]==1]
  }
  
  Y.mu <- apply(Y.raw,2,mean)
  Y.sd <- apply(Y.raw,2,sd)
  
  #browser()
  Y.raw <- apply(Y.raw, 2, function(x) (x-mean(x))/sd(x))
  
  X <- embed(Y.raw, dimension=p+1)
  X <- X[,(ncol(Y.raw)+1):ncol(X)]
  Y <- Y.raw[(p+1):nrow(Y.raw),]
  mid <- mid[(p+1):nrow(mid),]
  K <- ncol(X)
  
  # Some starting values for a linear VAR
  A.OLS <- solve(crossprod(X))%*%crossprod(X, Y)
  Sigma.OLS <- crossprod(Y-X%*%A.OLS)/nrow(X)
  
  Sigt <- array(0,dim=c(M,M,N))
  Sigt[,,1:N] <- Sigma.OLS #CHG
  Sig_big <- array(0,dim=c(M*p,M*p,N))
  Sig_big[1:M,1:M,1:N] <- Sigma.OLS #CHG
  
  A.prior <- A.OLS*0
  theta <- A.OLS^0 # shrinkage prior on the linear VAR coefficients
  a0sig <- prior.sig[1]
  b0sig <- prior.sig[2]
  
  # state eq priors
  beta0 <- Y[p,]
  for(j in 1:(p-1)){
    beta0 <- c(beta0,Y[p-j,])
  }
  
  pr_init <- 1
  P00 <- diag(length(beta0))*pr_init
  
  # Stochastic volatility
  priormu <- c(0, 1e-25)
  priorphi <- c(25, 5)
  priorsigma <- 1e-6
  
  # Initialize HS prior on covariances (if any)
  lambda.A <- 1
  nu.A <- 1
  tau.A <- 1
  zeta.A <- 1
  prior.cov <- rep(1, M*(M-1)/2)
  
  # initialize HS prior on linear VAR coefs (if VAR.mean=="linear")
  lambda.VAR <- 1
  nu.VAR <- 1
  tau.VAR <- 1
  zeta.VAR <- 1
  
  if(VAR.mean=="bart"){
    control <- dbartsControl(verbose = FALSE, keepTrainingFits = TRUE, useQuantiles = FALSE,
                             keepTrees = TRUE, n.samples = ntot,
                             n.cuts = 100L, n.burn = nburn, n.trees = num.trees, n.chains = 1,
                             n.threads = 1, n.thin = 1L, printEvery = 1,
                             printCutoffs = 0L, rngKind = "default", rngNormalKind = "default",
                             updateState = FALSE)
    sampler.list <- list()
    svdraw.list <- list()
    for (jj in seq_len(M)){
      sampler.list[[jj]] <- dbarts(Y[,jj]~X, control = control,tree.prior = cgm(cgm.exp, cgm.level), node.prior = normal(sd.mu),n.samples = nsave, weights=rep(1,N), sigma=sqrt(Sigma.OLS[jj,jj]), resid.prior = chisq(prior.sig[[1]], prior.sig[[2]]))
      svdraw.list[[jj]]  <- list(para = c(mu = -10, phi = 0.9, sigma = 0.2), latent = rep(0, N))
    }
    sampler.run <- list()  
  }
  
  Y0 <- Y
  Y0[mid[,1]==1,1] <- NA
  
  # Store some objects
  varcount <- matrix(NA_integer_, nsave, ncol(X))
  h.store <- array(NA, c(nsave, N, M))
  sigma2.store <- array(NA, c(nsave, M))
  A0.store <- array(NA, c(nsave, M, M))
  A.store <- array(NA, c(nsave, M, M*p))
  f.store <- array(NA, c(nsave, N, M))
  count.store <- array(NA, c(nsave, M, M*p))
  Ym.store <- array(NA,c(nsave,N))
  if(p==3){
    Yq.store <- array(NA,c(nsave,N/3)) 
  } else if(p==5){
    Yq.store <- array(NA,c(nsave,(N+2)/3))
  }

  #Some matrices to store estimates across equations
  sigma.mat <- matrix(NA, M, 1)
  A0.draw <- diag(M)
  eta.shocks <- matrix(NA, N, M)
  if(VAR.mean=="bart"){
    H.mat <- matrix(NA, N, M)  
  }else if(VAR.mean=="linear"){
    H.mat <- log(matrix(1, N, M))
  }
  HS.mat <- matrix(1, M, M)
  A <- matrix(0, M*p, M)
  count.mat <- matrix(0, M*p, M)
  
  # show progress
  if(!quiet){
    pb <- txtProgressBar(min = 0, max = ntot, style = 3)
    start  <- Sys.time()
  }
  
  # start sampling loop
  for (irep in seq_len(ntot)){
    count.var <- 0
    var.check <- TRUE
    while(var.check){
      count.var <- count.var + 1 
      if(VAR.mean=="bart"){
        # Simulate the constant error variance and the BART part
        X.ginv <- MASS::ginv(X)

        for (nr in seq_len(M)){
          if (nr > 1){
            Z <- eta.shocks[,1:(nr-1), drop=F]
            A0.nr <- A0.draw[nr,1:(nr-1)]
            sampler.list[[nr]]$setResponse(Y[,nr] - Z%*%A0.nr)
          }
          rep.i <- sampler.list[[nr]]$run(0L, 1L)
          
          sampler.run[[nr]] <- rep.i
          sigma.mat[nr,] <- rep.i$sigma
          if (any(is.na(rep.i$train))) break
          eta.shocks[,nr] <- Y[,nr] - rep.i$train 
          A[ , nr] <- X.ginv%*%rep.i$train
          count.mat[, nr] <- rep.i$varcount
          if (nr > 1){
            #Sample covariances 
            norm.nr <- as.numeric(exp(-.5*latent(svdraw.list[[nr]])) * 1/sigma.mat[nr,])
            u.nr <- eta.shocks[,1:(nr - 1), drop = FALSE] * norm.nr
            eta.nr <- eta.shocks[,nr] * norm.nr
            
            if (nr == 2) v0.inv <- 1/HS.mat[nr, 1] else v0.inv <- diag(1/HS.mat[nr, 1:(nr-1)])
            
            V.cov <- solve(crossprod(u.nr) + v0.inv)
            mu.cov <- V.cov %*% crossprod(u.nr, eta.nr)
            mu.cov.draw <- mu.cov + t(chol(V.cov)) %*% rnorm(ncol(V.cov)) 
            A0.draw[nr, 1:(nr-1)] <- mu.cov.draw
          }
          if (nr==1) gdp.fit <- (rep.i$train + rep.i$sigma*rnorm(N))*Y.sd[nr] + Y.mu[nr]
          if (exact & irep > nburn){
            f.store[irep-nburn,,nr] <- (rep.i$train + rep.i$sigma*rnorm(N))*Y.sd[nr] + Y.mu[nr]
          } 
        }

      }else if(VAR.mean=="linear"){
        for (nr in seq_len(M)){
          # mean coefficients
          norm <- exp(-H.mat[,nr]/2)
          XX <- X*norm
          if(nr == 1){
            YY <- Y[,nr]*norm
          }else{
            Z <- eta.shocks[,1:(nr-1), drop=F]
            A0.nr <- A0.draw[nr,1:(nr-1)]
            YY <- (Y[,nr] - Z%*%A0.nr)*norm
          }
          
          V_post <- try(solve((crossprod(XX)) + diag(1/theta[,nr])),silent=T)
          if (is(V_post,"try-error")) V_post <- ginv((crossprod(XX) + diag(1/theta[,nr])))
          A_mean <- V_post %*% (crossprod(XX, YY)+diag(1/theta[,nr])%*%A.prior[,nr])
          A.draw.nr <- try(A_mean+t(chol(V_post))%*%rnorm(K,0,1),silent=T)
          if (is(A.draw.nr,"try-error")) A.draw.nr <- mvrnorm(1,A_mean,V_post)
          A[,nr] <- A.draw.nr
          
          # covariances
          eta.shocks[,nr] <- Y[,nr] - X%*%A.draw.nr 
          if(nr>1){
            u.nr <- eta.shocks[,1:(nr-1),drop=F]*norm
            eta.nr <- eta.shocks[,nr]*norm
            
            if(nr==2){
              v0.inv <- 1/HS.mat[nr, 1]
            }else{
              v0.inv <- diag(1/HS.mat[nr, 1:(nr-1)])
            }
            
            V.cov <- solve(crossprod(u.nr) + v0.inv)
            mu.cov <- V.cov %*% crossprod(u.nr, eta.nr)
            mu.cov.draw <- mu.cov + t(chol(V.cov)) %*% rnorm(ncol(V.cov)) 
            A0.draw[nr, 1:(nr-1)] <- mu.cov.draw
          }
        }
      }
      
      # Simulate the stochastic volatility part conditional on the fit from the tree part of the model
      shock.normalized <- eta.shocks %*% t(solve(A0.draw))
      if(VAR.mean=="bart"){
        for (nr in seq_len(M)){
          svdraw.list[[nr]] <- svsample2(shock.normalized[ , nr] / sigma.mat[nr], startpara = para(svdraw.list[[nr]]), startlatent = latent(svdraw.list[[nr]]), priormu = priormu,  priorphi = priorphi, priorsigma = priorsigma)
          normalizer <- as.numeric(exp(-.5*latent(svdraw.list[[nr]])))
          weights.new <- as.numeric(exp(-latent(svdraw.list[[nr]])))
          dat <- dbartsData(formula = Y[, nr]~X,weights=weights.new)
          sampler.list[[nr]]$setData(dat)
          H.mat[, nr] <- log(sigma.mat[nr]^2) + latent(svdraw.list[[nr]]) 
        }
      }else if(VAR.mean=="linear"){
        for (nr in seq_len(M)){
          H.mat[,nr] <- log(1/rgamma(1,a0sig+N/2,b0sig+crossprod(shock.normalized[,nr])/2))
        }
      }
      
      # Sample scaling parameters of the Horseshoe prior for the full VC matrix
      hs_draw <- get.hs(bdraw=A0.draw[lower.tri(A0.draw)],lambda.hs=lambda.A,nu.hs=nu.A,tau.hs=tau.A,zeta.hs=zeta.A)
      lambda.A <- hs_draw$lambda
      nu.A <- hs_draw$nu
      tau.A <- hs_draw$tau
      zeta.A <- hs_draw$zeta
      prior.cov <- hs_draw$psi
      HS.mat[lower.tri(HS.mat)] <- prior.cov
      
      if(VAR.mean=="linear"){
        hs_draw <- get.hs(bdraw=as.numeric(A),lambda.hs=lambda.VAR,nu.hs=nu.VAR,tau.hs=tau.VAR,zeta.hs=zeta.VAR)
        lambda.VAR <- hs_draw$lambda
        nu.VAR <- hs_draw$nu
        tau.VAR <- hs_draw$tau
        zeta.VAR <- hs_draw$zeta
        theta <- matrix(hs_draw$psi,K,M)
      }
      
      # Interpolate missing values using FFBS
      MM <- get_companion(A,c(M,0,p))$MM
      
      for(tt in seq_len(N)){
        S_tmp <- exp(H.mat[tt,])
        S.t <- t(A0.draw)%*%crossprod(diag(S_tmp),(A0.draw))

        Sigt[,,tt] <- S.t
        Sig_big[1:M,1:M,tt] <- Sigt[,,tt]
      }
      
      KF.smooth.prop <-fkf(beta0, P00, matrix(0, K, 1), matrix(0, M, 1), MM, Ft, Sig_big, Rt, t(Y0), check.input = TRUE)
      if (any(is.na(KF.smooth.prop$att))){
        Y0.zero <- Y0
        Y0.zero[is.na(Y0.zero)] <- 0
        KF.smooth <- KF.slow(y=Y0.zero,MM = MM, Ft = Ft, Rt=Rt, Sig_big = Sig_big, kk = K, t = N, beta0 = beta0, P00 = P00)
      }else{
        KF.smooth <- KF.smooth.prop
      }
      
      # smoothing (backward recursions)
      beta2 <- t(KF.smooth$att)
      bm2 <- beta2
      jv <- 1:M # companion selector
      jv1 <- seq(1,K,by=M) # selector of lagged states
      
      ptt <- aperm(KF.smooth$Ptt, c(3,1,2))
      beta_tt <- aperm(KF.smooth$att, c(2,1))
      
      p00 <- ptt[N,jv1,jv1]
      beta2[N,] <- beta_tt[N,]
      beta2[N,jv1] <-  t(mvtnorm::rmvnorm(1, beta_tt[N, jv1], as.matrix(Matrix::forceSymmetric(p00))))
      
      q <- Sig_big[jv,jv,N]
      f <- MM[jv,]
      
      ztt <- matrix(0, N+p, 1)
      for (i in (N-1):1) {
        pt <- ptt[i,,]
        fpfq_inv <- try(solve(f%*%tcrossprod(pt,f)+q), silent=TRUE)
        if (is(fpfq_inv, "try-error")) fpfq_inv <- MASS::ginv(f%*%tcrossprod(pt,f)+q)
        
        bm <- beta_tt[i,]+t(pt%*%t(f)%*%fpfq_inv%*%t(beta2[i+1,jv]-beta_tt[i,]%*%t(f)))
        pm <- pt-pt%*%t(f)%*%fpfq_inv%*%f%*%pt
        beta2[i,] <- bm
        beta2.draw <- try(bm[jv1]+t(chol(pm[jv1,jv1]))%*%rnorm(p), silent=TRUE)
        if (is(beta2.draw, "try-error")) beta2.draw <- t(mvtnorm::rmvnorm(1, bm[jv1], as.matrix(Matrix::forceSymmetric(pm[jv1, jv1]))))
        beta2[i,jv1] <- beta2.draw
        bm2[i,] <- bm
      }
      
      if(p==3){
        ztt1 <- NULL
        for (t in seq(p, N+p)){
          ztt1 <- c(ztt1, beta2[t-p, seq(1, K, M)])
        }
        ztt1 <- ztt1[seq(1, N*p, p)]

        ztt <- NULL
        z.quarterly <- NULL
        for (t in seq(p, N+p, by=3)){
          ztt <- c(ztt, beta2[t-p,seq(1, K, M)])
          z.quarterly <- c(z.quarterly, mean(beta2[t-p,seq(1, K, M)][1:3]))
        }
        z.quarterly <- z.quarterly[-1]
      }else if(p==5){
        sl.agg <- mid.short[,1]==0; sl.agg[length(sl.agg)] <- TRUE
        zt.mat <- beta2[sl.agg,seq(1,K,M)]
        ztt <- as.vector(t(zt.mat[,1:3]))

        #Compute difference between ztt and N
        N.diff <- length(ztt) - N
        ztt <- ztt[(N.diff+1):length(ztt)]
        
        z.quarterly <- apply(loadings*t(zt.mat),2,sum)
      }
      
      var.check <- !(sd(beta2[,1]) < var.thrsh)
      if (count.var == max.count.var){
        ztt <- (X%*%A)[,1]
        message("Reached maximum amount of replications")
        count.var <- FALSE
        var.check <- FALSE
      }
    }
    
    X.mm <- Y
    X.mm[1:N, ] <- cbind(ztt,Y[,2:M])
    out.prop <- rbind(matrix(0, p, M), X.mm)
    out <- out.prop
    
    Y <- out
    X <- embed(Y, dimension=p+1)
    X <- X[, (ncol(Y)+1):ncol(X)]
    Y <- Y[(p+1):nrow(Y),]
    
    # Final step: Store draws of interest
    if (irep > nburn & sd(ztt)<var.thrsh){
      sigma2.store[irep-nburn,] <- t(sigma.mat^2)
      A.store[irep-nburn,,] <- t(A)
      A0.store[irep-nburn,,] <- A0.draw
      h.store[irep-nburn,,] <- H.mat
      if(exact){
        Ym.store[irep-nburn,] <- f.store[irep-nburn,,]
      }else{
        Ym.store[irep-nburn,] <- (Y[,1]*Y.sd[1])+Y.mu[1]
      }
      Yq.store[irep-nburn,] <- (z.quarterly*Y.sd[1])+Y.mu[1]
      count.store[irep-nburn,,] <- t(count.mat)
    }
    if(!quiet){
      setTxtProgressBar(pb, irep)
      if (irep %% iter.update==0){
        end <- Sys.time()
        message(paste0("\n Average time necessary to produce a single draw over the last ",iter.update," draws ", round(as.numeric(end-start)/iter.update, digits=4), " seconds, currently at draw ", irep,"."))
        start <- Sys.time()
        par(mfrow=c(1,1))
        ts.plot(Y[,1])
        lines(gdp.fit)
        points(Y0[,1])
        
      }
    }
  }
  
  return(list("A"=A.store,"A0"=A0.store,"sigma2"=sigma2.store,"h"=h.store,
              "data_m"=Ym.store, "fitvals"=f.store, "X"=X))
}
