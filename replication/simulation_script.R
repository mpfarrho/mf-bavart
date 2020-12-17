require(dbarts)
require(stochvol)
require(RColorBrewer)
require(snowfall)
require(parallel)
require(bvarsv)
require(mfbvar)
require(FKF)
require(Matrix)

source("mfbart.R")

set.seed(123)
#Simulate a mixed frequency VAR model
replications <- 150
intensity <- c(0, 6, 15, 18, 25)
start.crisis <- c(10,20, 25, 35)
dgp <- c("VAR", "NL")
expand.i <- expand.grid(intensity, start.crisis, dgp, 1:replications)
run <- 1

sl.i <- expand.i[run, ]

intensity <- 6#sl.i[[1]]
length.crisis <- 10#sl.i[[2]]
dgp <- as.character(sl.i[[3]])
scale.crisis <- c(0, intensity)
rep.i <- sl.i[[4]]

dgp.set <- "VAR"
T <- 251
M <- 5
p.true <- 5
Y.latent <- matrix(0, T, M)
Y.observed <- matrix(NA, T, M)
n.quarter <- 1
#Start by simulating the latent process (which is a standard VAR model)
if (dgp.set == "VAR"){
  stab.cond <- TRUE
  while (stab.cond){
    A.true <- matrix(rnorm(M*M, 0, 0.1), M, M)
    diag(A.true) <- 0.9
    
    if (max(abs(Re(eigen(A.true)$values))) < 1) stab.cond <- FALSE
    
    Q.true <- diag(M)*0.4
    Q.true[lower.tri(Q.true)] <- rnorm(M*(M-1)/2, 0, 0.01)
    Sigma.true <- crossprod(Q.true)
    Y.latent[1,] <- rnorm(M,0, 1)
    for (t in 2:T){
      Y.latent[t,] <- A.true %*% Y.latent[t-1,] + Q.true%*%rnorm(M)
    }
  }
  sd.Y <- sd(Y.latent[,1])
  
  Finv <- solve(diag(M) - A.true)
  Sig.un <- sqrt(diag(Finv %*% Sigma.true %*% t(Finv)))
  
   scale.recovery <- c(2, 3)
   if (intensity > 0){
    end.crisis <- 4
  
     Y.latent[((T-length.crisis+1):(T)), 1] <- -sd.Y*seq(scale.crisis[[1]], scale.crisis[[2]], length.out = length.crisis)
     Y.latent[(T-2+1):(T), 1]   <- sd.Y*seq(scale.recovery[[1]], scale.recovery[[2]], length.out = 2)
   }
   start.crisis <- length.crisis
}else{
  #Non-linear function
  A.true <- matrix(rnorm(M*M, 0, 0.1), M, M)
  diag(A.true) <- 0

  Q.true <- diag(M)*1
  Q.true[lower.tri(Q.true)] <- rnorm(M*(M-1)/2, 0, 0.01)
  Sigma.true <- crossprod(Q.true)
  Y.latent[1:2,] <- rnorm(2*M,2, 1)
  
  for (t in 3:T){
    Y.latent[t,] <-  A.true %*% Y.latent[t-1, ]/as.numeric(Y.latent[t-2,3]) + matrix(rnorm(M*M, 0, 1e-2), M, M)%*%(Y.latent[t-1,] * as.numeric(Y.latent[t-1, 4])) + Q.true%*%rnorm(M)#(A.true %*% Y.latent[t-1,])*ind.i + (A.true.2 %*% Y.latent[t-1,]) *(1-ind.i) + Q.true%*%rnorm(M)
  }
  sd.Y <- sd(Y.latent[,1])
  end.crisis <- 0
  ts.plot(Y.latent[,1])
}

sl.quarter <- seq(5,T, by=3)
for (i in 1:T){
  if (i %in% sl.quarter){
    load <- 1/9*c(1,rep(0, M-1), 2, rep(0, M-1), 3, rep(0, M-1), 2, rep(0, M-1), 1,rep(0, M-1))
    Y.lags <- c(Y.latent[i,],Y.latent[i-1,], Y.latent[i-2,], Y.latent[i-3,], Y.latent[i-4,])
    Y.obs <- load%*%Y.lags
    Y.observed[i,1] <- Y.obs
  }
}
Y.observed[,2:M] <- Y.latent[, 2:M]

nsave <- 500 #Save 1000 draws
nburn <- 500 # Burn 1000 draws
p <- 5 # laglength
prior.cov <- 0.01
sd.mu <- 2
var.thrsh <- 3
count.var.max <- 100

cgm.level <- 0.95
cgm.exp <- 2
num.trees <- 250

ntot <- nsave+nburn
prior.sig <- c(nrow(Y.observed)/2, 0.75)
mf.country <- mfbart(Y.observed,p=p,VAR.mean = "bart",nburn=nburn,nsave=nsave,prior.cov=0.01,cgm.level=cgm.level,cgm.exp=cgm.exp,sd.mu=sd.mu,num.trees=num.trees,prior.sig = prior.sig,n.quarter=1,quiet=FALSE,iter.update=50,var.thrsh=3, max.count.var = 20, exact=FALSE)
mf.var <- mfbart(Y.observed,p=p,VAR.mean = "linear",nburn=nburn,nsave=nsave,prior.cov=0.01,cgm.level=cgm.level,cgm.exp=cgm.exp,sd.mu=sd.mu,num.trees=num.trees,prior.sig = prior.sig,n.quarter=1,quiet=FALSE,iter.update=50,var.thrsh=3, max.count.var = 20, exact=FALSE)

pred.dens.BART <- t(apply(mf.country$data_m,2,function(x) quantile(x, c(0.05 ,0.5,0.95),na.rm=T)))
pred.dens.VAR <- t(apply(mf.var$data_m,2,function(x) quantile(x, c(0.16 ,0.5,0.84), na.rm=T)))

foldername <- "Results"
dir.create(foldername, showWarnings = FALSE)

pdf(paste0(foldername,"/one_sim", dgp.set, ".pdf"), width=16, height=8)
Y.out <- Y.latent[5:(T-1)]
ts.plot(Y.out, col=c(0,1,1,2),  ylab ="", xlab="")
polygon(c(1:nrow(pred.dens.BART),rev(1:nrow(pred.dens.BART))),c(pred.dens.BART[,1],rev(pred.dens.BART[,3])),col="lightgray",border=NA)##,irf_true[jj,kk,]
lines(pred.dens.BART[,1]);lines(pred.dens.BART[,3]);lines(pred.dens.BART[,2])
lines(Y.out, col = "red")
dev.off()
