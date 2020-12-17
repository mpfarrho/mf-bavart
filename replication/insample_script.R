source("mfbart.R")
nsave <- 5000 #Save 1000 draws
nburn <- 5000 # Burn 1000 draws

p <- 5 # laglength

# stability checks
sd.mu <- 2
var.thrsh <- 3
count.var.max <- 100

cgm.level <- 0.95 # alpha
cgm.exp <- 2 # beta
num.trees <- 250 # S
prior.cov <- 0.01

ntot <- nsave+nburn

data.all <- FALSE # use pandemic dataset
cn.run <- "DE" # choose from c("IT","DE", "FR", "ES")
source("data_script.R")
prior.sig <- c(nrow(Y.m)/2, 0.75)
mf.country <- mfbart(Y.m,p=p,VAR.mean="bart",nburn=nburn,nsave=nsave,prior.cov=0.01,cgm.level=cgm.level,cgm.exp=cgm.exp,sd.mu=sd.mu,num.trees=num.trees,prior.sig = prior.sig,n.quarter=1,quiet=FALSE,iter.update=50,var.thrsh=3,max.count.var=100,exact=TRUE)
