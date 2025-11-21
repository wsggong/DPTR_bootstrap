FD_transform <- function(y, n, t){
  FDy <- matrix(NA, nrow=n, ncol=t)
  Y <- matrix(y, nrow=n, ncol=t)
  for(i in 2:t){
    FDy[,i] <- Y[,i]-Y[,i-1]
  }
  return(c(FDy))
}

lag_dpd <- function(y, n, t, s){
  
  res <- rep(NA, n*t)
  res[(n*s+1):(n*t)] <- y[1:(n*(t-s))]
  
  return(res)
}

lag_dpd_inst <- function(y, n, t, lb=NULL, ub=NULL, t0=NULL, cut=FALSE, gamma=NULL, q=NULL){
  
  if(is.null(ub)){
    ub <- (t-1)
  }else{
    ub <- min(t-1, ub)
  }
  
  if(is.null(t0)){
    t0 <- lb+1
  }
  
  if(t0-lb <=0){
    stop("Error, t0-lb > 0 should hold")
  }
  
  if(lb==ub){
    dimz <- t-t0+1
  }else{
    dimz <- (ub-t0+1)*(ub+t0)/2-lb*(ub-t0+1)+(t-ub)*(ub-lb+1)
  }
  
  res <- matrix(0, nrow=n*(t-t0+1), ncol=dimz)
  Y <- matrix(y, nrow=n, ncol=t)
  iold <- inew <- 1
  
  for(i in t0:t){
    res_t <- Y[,max(i-ub,1):(i-lb)]
    inew <- iold+((i-lb)-max(i-ub,1))
    res[(n*(i-t0)+1):(n*(i-t0+1)),iold:inew] <- res_t
    iold <- inew+1
  }
  nix <- !is.na(colSums(res))
  return(res[,nix])
  
}

arrange_inst <- function(Z, n, t, t0, collapse=FALSE){
  ix <- matrix(FALSE, nrow=(t-t0+1), ncol=ncol(Z))
  
  for(i in 1:(t-t0+1)){
    for(j in 1:ncol(Z)){
      ix[i,j] <- sum(Z[(n*(i-1)+1):(n*i),j]!=0)>0
    }
  }
  
  if(collapse==TRUE){
    res <- matrix(0, nrow=nrow(Z), ncol=max(rowSums(ix)))
    for(i in 1:(t-t0+1)){
      res[(n*(i-1)+1):(n*i),(ncol(res)-rowSums(ix)[i]+1):ncol(res)] <- Z[(n*(i-1)+1):(n*i),ix[i,]]
    }
  }else{
    res <- matrix(0, nrow=nrow(Z), ncol=ncol(Z))
    nix <- c(0,cumsum(rowSums(ix)))
    for(i in 1:(t-t0+1)){
      res[(n*(i-1)+1):(n*i),(nix[i]+1):nix[i+1]] <- Z[(n*(i-1)+1):(n*i),ix[i,]]
    }
  }
  
  return(res)
}



#========================================#
#        Unconstrained Estimation        #
#========================================#

gmm_fixed_est <- function(y,                     # Dependent Variable
                          x1,                    # Regressor 1
                          x2,                    # Regressor 2
                          q,                     # Threshold variable
                          gamma,                 # Fixed Threshold Location
                          weight=NULL,           # GMM weight matrix
                          weight_out=FALSE,      # Compute Omega inverse? T/F
                          inst,                  # Customized Instrument
                          n=NULL,                # Cross sectional sample size
                          t=NULL,                # Time length
                          t0=NULL                # First time used
                          ){
  
  if((is.null(n)|is.null(t))){
    stop("The dimension of panel data should be determined either by providing 'n' and 't'")
  }
  if(is.null(inst)){
    stop("'inst' is not given.")
  }
  
  p1 <- length(x1)/(n*t)
  p2 <- length(x2)/(n*t)
  
  x1 <- as.matrix(x1)
  x2 <- as.matrix(x2)
  FDy <- FD_transform(y, n, t)
  FDx <- matrix(NA, nrow=n*t, ncol=p1)
  for(i in 1:p1){
    FDx[,i] <-FD_transform(x1[,i], n, t)
  }
  
  Y <- FDy[(n*(t0-1)+1):(n*t)]
  
  X1 <- FDx[(n*(t0-1)+1):(n*t),]
  X2 <- cbind(1*(q[(n*(t0-1)+1):(n*t)]>gamma)-1*(q[(n*(t0-2)+1):(n*(t-1))]>gamma),
              x2[(n*(t0-1)+1):(n*t),]*(q[(n*(t0-1)+1):(n*t)]>gamma)-x2[(n*(t0-2)+1):(n*(t-1)),]*(q[(n*(t0-2)+1):(n*(t-1))]>gamma))
    
  
  X <- cbind(X1,X2)
  Z <- inst
  
  if(is.null(weight)){
    W <- diag(rep(1,ncol(Z)))
    # Z_a <- rbind(Z, matrix(0,nrow=n,ncol=ncol(Z)))-rbind(matrix(0,nrow=n,ncol=ncol(Z)),Z)
    # W <- MASS::ginv(t(Z_a)%*%Z_a)/n
  }else{
    W <- weight
  }
  
  g1_bar <- t(Z) %*% (Y)/n
  g2_bar <- t(Z) %*% (X)/n
  
  singular <- FALSE
  tryCatch(
    coef <- solve(t(g2_bar) %*% W %*% g2_bar, t(g2_bar) %*% W %*% g1_bar),
    error=function(e){singular <<- TRUE}
  )
  
  if(singular==FALSE){
    g_bar <- g1_bar - g2_bar %*% coef          # GMM objective
    obj <- t(g_bar) %*% W %*% g_bar   # GMM objective
    
    if(weight_out==TRUE){
      ghat_full <- Z * c(Y-X%*%coef)
      ghat <- matrix(0, nrow=n, ncol=ncol(Z))
      
      for(i in 1:(t-(t0-1))){
        ghat <- ghat + ghat_full[(n*(i-1)+1):(n*i),]
      }; rm(ghat_full)
      W <- MASS::ginv(cov(ghat))
    }
    
    return(list(coefficients=coef, objective=obj, weight=W, g_bar=g_bar))
  }else{
    return(list(coefficients=NA, objective=+Inf, weight=NA, g_bar=NA))
  }
}



gmm_unconstrained_est <- function(y, 
                                  x1, 
                                  x2, 
                                  q, 
                                  q_grid, 
                                  two_step=TRUE, 
                                  weight=NULL, 
                                  weight_out=FALSE, 
                                  instmethod=1, 
                                  inst=NULL,
                                  n=NULL, 
                                  t=NULL,
                                  t0=3){
  
  objlist <- rep(NA, length(q_grid))
  
  # First step
  for(q_iter in 1:length(q_grid)){
    gmm0 <- gmm_fixed_est(y, x1, x2, q, q_grid[q_iter], weight, weight_out=F, inst, n, t, t0)
    objlist[q_iter] <- gmm0$objective
  }
  
  # Second Step
  if(two_step==TRUE){
    ix <- which.min(objlist)
    weight <- gmm_fixed_est(y, x1, x2, q, q_grid[ix], weight, weight_out=T, inst, n, t, t0)$weight
    for(q_iter in 1:length(q_grid)){
      gmm0 <- gmm_fixed_est(y, x1, x2, q, q_grid[q_iter], weight, weight_out=F, inst, n, t, t0)
      objlist[q_iter] <- gmm0$objective
    }
  }
  
  ix <- which.min(objlist)
  gmm <- gmm_fixed_est(y, x1, x2, q, q_grid[ix], weight, weight_out=weight_out, inst, n, t, t0)
  
  return(list(coefficients=gmm$coefficients, 
              threshold=q_grid[ix], 
              objective=gmm$objective, 
              objlist=objlist, 
              used_grid=q_grid, 
              weight=gmm$weight,
              g_bar=gmm$g_bar))
}




#========================================#
#    Continuity Restricted Estimation    #
#========================================#
# 
gmm_fixed_constrained_est <- function(y,                     # Dependent Variable
                                      x1,                    # Regressor 1
                                      q,                     # Threshold variable
                                      gamma,                 # Fixed Threshold Location
                                      weight=NULL,           # GMM weight matrix
                                      weight_out=FALSE,      # Compute Omega inverse? T/F
                                      inst,                  # Customized Instrument
                                      n=NULL,                # Cross sectional sample size
                                      t=NULL,                # Time length
                                      t0=NULL                # First time used
                                      ){
  if((is.null(n)|is.null(t))){
    stop("The dimension of panel data should be determined either by providing 'n' and 't'")
  }
  if(is.null(inst)){
    stop("'inst' is not given.")
  }
  
  p1 <- length(x1)/(n*t)
  
  x1 <- as.matrix(x1)
  FDy <- FD_transform(y, n, t)
  FDx <- matrix(NA, nrow=n*t, ncol=p1)
  for(i in 1:p1){
    FDx[,i] <-FD_transform(x1[,i], n, t)
  }
  
  Y <- FDy[(n*(t0-1)+1):(n*t)]
  
  X1 <- FDx[(n*(t0-1)+1):(n*t),]
  X2 <- (q[(n*(t0-1)+1):(n*t)]-gamma)*(q[(n*(t0-1)+1):(n*t)]>gamma)-(q[(n*(t0-2)+1):(n*(t-1))]-gamma)*(q[(n*(t0-2)+1):(n*(t-1))]>gamma)
  
  X <- cbind(X1,X2)
  Z <- inst
  
  if(is.null(weight)){
    W <- diag(rep(1,ncol(Z)))
    # Z_a <- rbind(Z, matrix(0,nrow=n,ncol=ncol(Z)))-rbind(matrix(0,nrow=n,ncol=ncol(Z)),Z)
    # W <- MASS::ginv(t(Z_a)%*%Z_a)/n
  }else{
    W <- weight
  }
  
  g1_bar <- t(Z) %*% (Y)/n
  g2_bar <- t(Z) %*% (X)/n
  
  singular <- FALSE
  tryCatch(
    coef <- solve(t(g2_bar) %*% W %*% g2_bar, t(g2_bar) %*% W %*% g1_bar),
    error=function(e){singular <<- TRUE}
  )
  
  if(singular==FALSE){
    g_bar <- g1_bar - g2_bar %*% coef          # GMM objective
    obj <- t(g_bar) %*% W %*% g_bar   # GMM objective
    
    if(weight_out==TRUE){
      ghat_full <- Z * c(Y-X%*%coef)
      ghat <- matrix(0, nrow=n, ncol=ncol(Z))
      
      for(i in 1:(t-(t0-1))){
        ghat <- ghat + ghat_full[(n*(i-1)+1):(n*i),]
      }; rm(ghat_full)
      W <- MASS::ginv(cov(ghat))
    }
    return(list(coefficients=coef, objective=obj, weight=W, g_bar=g_bar))
  }else{
    return(list(coefficients=NA, objective=+Inf, weight=NA, g_bar=NA))
  }
}


gmm_constrained_est <-  function(y, 
                                 x1, 
                                 q, 
                                 q_grid, 
                                 two_step=TRUE, 
                                 weight=NULL, 
                                 weight_out=FALSE, 
                                 instmethod=1, 
                                 inst=NULL,
                                 n=NULL, 
                                 t=NULL,
                                 t0=3){
  
  objlist <- rep(NA, length(q_grid))
  
  
  # First step
  for(q_iter in 1:length(q_grid)){
    gmm0 <- gmm_fixed_constrained_est(y, x1, q, q_grid[q_iter], weight, weight_out=F, inst, n, t, t0)
    objlist[q_iter] <- gmm0$objective
  }
  
  # Second Step
  if(two_step==TRUE){
    ix <- which.min(objlist)
    weight <- gmm_fixed_constrained_est(y, x1, q, q_grid[ix], weight, weight_out=T, inst, n, t, t0)$weight
    for(q_iter in 1:length(q_grid)){
      gmm0 <- gmm_fixed_constrained_est(y, x1, q, q_grid[q_iter], weight, weight_out=F, inst, n, t, t0)
      objlist[q_iter] <- gmm0$objective
    }
  }
  ix <- which.min(objlist)
  gmm <- gmm_fixed_constrained_est(y, x1, q, q_grid[ix], weight, weight_out=weight_out, inst, n, t, t0)
  
  
  return(list(coefficients=gmm$coefficients, 
              threshold=q_grid[ix], 
              objective=gmm$objective, 
              objlist=objlist, 
              used_grid=q_grid, 
              weight=gmm$weight,
              g_bar=gmm$g_bar))
}



#==========================#
#    Design estimation     #
#==========================#

mNW <- function(x, X, Y, h=NULL, K = dnorm, stationary=FALSE) {
  
  # Arguments
  # x: evaluation points (gamma0)
  # X: vector (size n) with the predictors
  # Y: vector (size n) with the response variable
  # h: bandwidth
  # K: kernel
  
  # Silverman's rule of thumb
  if(is.null(h)){
    h <- bw.nrd(X)
  }
  n <- length(Y)
  
  # Matrix of size n x length(x)
  Kx <- sapply(X, function(Xi) K((x - Xi) / h)/h)
  
  if(stationary==TRUE){
    if(length(x)==1) W <- Kx / sum(Kx)
    else W <- Kx / rowSums(Kx)
    return(W %*% Y)
  }
  else{
    # Means at x ("drop" to drop the matrix attributes)
    return((Kx %*% Y)/n)
  }
  
}


Design_est <- function(Z, y, x1, x2, q, gamma, delta, delta3, n, t, t0=3, h=NULL){
  
  x1 <- as.matrix(x1)
  x2 <- as.matrix(x2)
  p1 <- length(x1)/(n*t)
  p2 <- length(x2)/(n*t)
  
  FDx <- matrix(NA, nrow=n*t, ncol=p1)
  for(i in 1:p1){
    FDx[,i] <-FD_transform(x1[,i], n, t)
  }
  
  X1 <- FDx[(n*(t0-1)+1):(n*t),]
  X2 <- cbind(1*(q[(n*(t0-1)+1):(n*t)]>gamma)-1*(q[(n*(t0-2)+1):(n*(t-1))]>gamma),
              x2[(n*(t0-1)+1):(n*t),]*(q[(n*(t0-1)+1):(n*t)]>gamma)-x2[(n*(t0-2)+1):(n*(t-1)),]*(q[(n*(t0-2)+1):(n*(t-1))]>gamma))
  
  X <- cbind(X1,X2)
  Mn <- -t(Z) %*% (X)/n
  ix <- matrix(FALSE, nrow=(t-t0+1), ncol=ncol(Z))
  for(i in 1:(t-t0+1)){
    for(j in 1:ncol(Z)){
      ix[i,j] <- sum(Z[(n*(i-1)+1):(n*i),j]!=0)>0
    }
  }
  nix <- c(0,cumsum(rowSums(ix)))
  Gn <- rep(NA, length=ncol(Z))
  for(i in 1:(t-t0+1)){
    for(j in (nix[i]+1):nix[i+1]){
      Gn[j] <- mNW(gamma, q[(n*(i+t0-2)+1):(n*(i+t0-1))], Z[(n*(i-1)+1):(n*i),j]*(cbind(1,x2[(n*(i+t0-2)+1):(n*(i+t0-1)),])%*%delta), h) - mNW(gamma, q[(n*(i+t0-3)+1):(n*(i+t0-2))], Z[(n*(i-1)+1):(n*(i)),j]*(cbind(1,x2[(n*(i+t0-3)+1):(n*(i+t0-2)),])%*%delta), h)
    }
  }
  
  Hn <- rep(NA, length=ncol(Z))
  for(i in 1:(t-t0+1)){
    for(j in (nix[i]+1):nix[i+1]){
      Hn[j] <- delta3/2*(mNW(gamma, q[(n*(i+t0-2)+1):(n*(i+t0-1))], Z[(n*(i-1)+1):(n*i),j], h) - mNW(gamma, q[(n*(i+t0-3)+1):(n*(i+t0-2))], Z[(n*(i-1)+1):(n*i),j], h))
    }
  }
  
  return(list(Mn=Mn, Gn=Gn, Hn=Hn))
}




gmm_linear_est <- function(y,                     # Dependent Variable
                          x1,                    # Regressor 1
                          weight=NULL,           # GMM weight matrix
                          inst,                  # Customized Instrument
                          n=NULL,                # Cross sectional sample size
                          t=NULL,                # Time length
                          t0=3                # First time used
){
  
  if((is.null(n)|is.null(t))){
    stop("The dimension of panel data should be determined either by providing 'n' and 't'")
  }
  if(is.null(inst)){
    stop("'inst' is not given.")
  }
  
  p1 <- length(x1)/(n*t)
 
  x1 <- as.matrix(x1)
  
  FDy <- FD_transform(y, n, t)
  FDx <- matrix(NA, nrow=n*t, ncol=p1)
  for(i in 1:p1){
    FDx[,i] <-FD_transform(x1[,i], n, t)
  }
  
  Y <- FDy[(n*(t0-1)+1):(n*t)]
  
  X <- FDx[(n*(t0-1)+1):(n*t),]
  
  Z <- inst
  
  if(is.null(weight)){
    W <- diag(rep(1,ncol(Z)))
    # Z_a <- rbind(Z, matrix(0,nrow=n,ncol=ncol(Z)))-rbind(matrix(0,nrow=n,ncol=ncol(Z)),Z)
    # W <- MASS::ginv(t(Z_a)%*%Z_a)/n
  }else{
    W <- weight
  }
  
  g1_bar <- t(Z) %*% (Y)/n
  g2_bar <- t(Z) %*% (X)/n
  
  coef <- solve(t(g2_bar) %*% W %*% g2_bar, t(g2_bar) %*% W %*% g1_bar)
  
  g_bar <- g1_bar - g2_bar %*% coef          # GMM objective
  obj <- t(g_bar) %*% W %*% g_bar   # GMM objective

  ghat_full <- Z * c(Y-X%*%coef)
  ghat <- matrix(0, nrow=n, ncol=ncol(Z))
  
  for(i in 1:(t-(t0-1))){
    ghat <- ghat + ghat_full[(n*(i-1)+1):(n*i),]
  }; rm(ghat_full)
  W <- MASS::ginv(cov(ghat))
  
  coef <- solve(t(g2_bar) %*% W %*% g2_bar, t(g2_bar) %*% W %*% g1_bar)
  
  g_bar <- g1_bar - g2_bar %*% coef          # GMM objective
  obj <- t(g_bar) %*% W %*% g_bar   # GMM objective
  
  return(list(coefficients=coef, objective=obj, weight=W, g_bar=g_bar))
}
