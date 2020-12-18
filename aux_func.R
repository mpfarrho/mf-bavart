# "Nowcasting in a Pandemic using Non-Parametric Mixed Frequency VARs" 
# Huber, F., Koop, G., Onorante, L., Pfarrhofer, M., and J. Schreiner, 
# Journal of Econometrics (forthcoming)
# 
# This file contains several auxiliary functions.

# get posteriors for the horseshoe prior (see Makalic & Schmidt, 2015)
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

# lag variables
mlag <- function(X,lag){
  p <- lag
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- ncol(X)
  Xlag <- matrix(NA,Traw,p*N)
  for (ii in 1:p){
    Xlag[(p+1):Traw,(N*(ii-1)+1):(N*ii)]=X[(p+1-ii):(Traw-ii),(1:N)]
  }
  return(Xlag)
}

# fill initial missings (from mfbvar, by S. Ankargren)
fill_na <- function(Y) {
  apply(Y, 2, function(x) {
    n_x <- length(x) # save lentgh
    if (any(is.na(x))) {
      x <- x[1:max(which(is.na(x) == FALSE))] # get rid of NAs in the end
      for (i in which(is.na(x))) {
        x1 <- NA
        counter <- 1
        while (is.na(x1) == TRUE) {
          x1 <- x[i + counter]
          counter <- counter + 1
        }
        x[i] <- x1
      }
      
      trimmed_length <- length(x)
      if (trimmed_length < n_x) {
        x <- c(x, rep(NA, n_x - trimmed_length))
        for (i in trimmed_length:n_x) {
          x[i] <- x[trimmed_length]
        }
      }
    }
    x})
}

# construct loading matrix (from mfbvar, by S. Ankargren)
build_Lambda2 <- function(aggregation, n_lags) {
  n_vars <- length(aggregation)
  if (any(aggregation %in% "grw") && n_lags < 5) {
    Lambda <- matrix(0, n_vars, n_vars * 5)
  } else if (n_lags > 2) {
    Lambda <- matrix(0, n_vars, n_vars * n_lags)
  } else {
    stop("Too few lags!")
  }
  
  n_pseudolags <- dim(Lambda)[2]/n_vars
  for (i in 1:n_vars) {
    if (aggregation[i] == "m") {
      fill_vec <- c(1, rep(0, n_pseudolags - 1))
    }
    if (aggregation[i] == "lvl") {
      fill_vec <- c(rep(1/3, 3), rep(0, n_pseudolags - 3))
    }
    if (aggregation[i] == "grw") {
      fill_vec <- c(1/3, 2/3, 1, 2/3, 1/3, rep(0, n_pseudolags - 5))/3 # Divide by three to make commensurate in scale
    }
    
    Lambda[i, seq(i, n_pseudolags * n_vars, by = n_vars)] <- fill_vec
  }
  return(Lambda)
}

# construct design matrix from ts-objects (from mfbvar, by S. Ankargren)
list_to_matrix <- function(Y_in) {
  require(lubridate)
  if (all(sapply(Y_in, function(x) inherits(x, "ts"))) || all(sapply(Y_in, function(x) inherits(x, "zoo")))) {
    if (all(sapply(Y_in, function(x) inherits(x, "ts")))) {
      zoofun <- function(x) {
        if (frequency(x) == 4) {
          if (is.null(dim(x))) {
            zoo::zoo(as.numeric(x), as.Date(zoo::as.Date.ts(x) %m+% months(2)))
          } else {
            zoo::zoo(as.matrix(x), as.Date(zoo::as.Date.ts(x) %m+% months(2)))
          }
        } else if (frequency(x) == 12) {
          if (is.null(dim(x))) {
            zoo::zoo(as.numeric(x), as.Date(zoo::as.Date.ts(x)))
          } else {
            zoo::zoo(as.matrix(x), as.Date(zoo::as.Date.ts(x)))
            
          }
        } else {
          stop("The data must only include monthly and/or quarterly time series.")
        }
      }
      
    } else if (all(sapply(Y_in, function(x) inherits(x, "zooreg")))) {
      zoofun <- function(x) {
        if (frequency(x) == 4) {
          if (is.null(dim(x))) {
            zoo::zoo(as.numeric(x), as.Date(zoo::as.Date(zoo::index(x)) %m+% months(2)))
          } else {
            zoo::zoo(as.matrix(x), as.Date(zoo::as.Date(zoo::index(x)) %m+% months(2)))
          }
        } else if (frequency(x) == 12) {
          if (is.null(dim(x))) {
            zoo::zoo(as.numeric(x), as.Date(zoo::as.Date(zoo::index(x))))
          } else {
            zoo::zoo(as.matrix(x), as.Date(zoo::as.Date(zoo::index(x))))
          }
        } else {
          stop("The data must only include monthly and/or quarterly time series.")
        }
      }
    }
    zoolist <- lapply(Y_in, zoofun)
    reducedlist <- Reduce(zoo::merge.zoo, zoolist)
    Y <- as.matrix(reducedlist)
    rownames(Y) <- as.character(time(reducedlist))
    dim_null <- sapply(zoolist, function(x) is.null(dim(x)))
    if (all(dim_null)) {
      colnames(Y) <- names(zoolist)
    } else if (all(!dim_null)) {
      colnames(Y) <- Reduce(c, lapply(zoolist, colnames))
    } else {
      name_vec <- c()
      for (iter in 1:length(dim_null)) {
        if (dim_null[iter]) {
          name_vec <- c(name_vec, names(zoolist)[iter])
        } else {
          name_vec <- c(name_vec, colnames(zoolist[[iter]]))
        }
      }
      colnames(Y) <- name_vec
    }
    
    if (all(dim_null)) {
      zoolistfreq <- sapply(Y_in, frequency)
    } else if (all(!dim_null)) {
      zoolistfreq <- sapply(Y_in, frequency)
      zoolistn <- sapply(Y_in, NCOL)
      zoolistfreq <- Reduce(c, mapply(function(x, y) rep(x, each = y), zoolistfreq, zoolistn, SIMPLIFY =  FALSE))
      
    } else {
      zoolistfreq <- c()
      for (iter in 1:length(dim_null)) {
        if (dim_null[iter]) {
          zoolistfreq <- c(zoolistfreq, frequency(Y_in[[iter]]))
        } else {
          zoolistfreq <- c(zoolistfreq, rep(frequency(Y_in[[iter]]), each = ncol(Y_in[[iter]])))
        }
      }
    }
    names(zoolistfreq) <- NULL
    if (all(zoolistfreq %in% c(4, 12))) {
      freq <- ifelse(zoolistfreq == 4, "l", "h")
    } else {
      stop("Only monthly and quarterly frequencies are allowed.")
    }
  } else {
    
  }
  return(list("Yraw"=Y, "freq"=freq))
}

# get VAR in companion form
get_companion <- function(Beta_,varndxv){
  nn <- varndxv[[1]]
  nd <- varndxv[[2]]
  nl <- varndxv[[3]]
  
  nkk <- nn*nl+nd
  
  Jm <- matrix(0,nkk,nn)
  Jm[1:nn,1:nn] <- diag(nn)
  
  MM <- rbind(t(Beta_),cbind(diag((nl-1)*nn), matrix(0, (nl-1)*nn, nn)))
  
  return(list(MM=MM,Jm=Jm))
}
