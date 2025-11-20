#=====================================#
#     Null imposed Bootstrapping      #
#=====================================#

boot_null <- function(Z, # IV
                      y,  # outcome
                      x1,  # regressor
                      x2,  # regressor, same with x1 in our case
                      q,  # threshold variable
                      alpha, #alpha0*
                      alpha_hat, #alpha_hat (for constructing residual)
                      gamma, #gamma0*
                      gamma_hat, #gamma_hat (for constructing residual)
                      n=NULL, 
                      t=NULL, 
                      t0=3,
                      rng=NULL){
  
  # 1. Obtain residuals
  
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
  X2 <- cbind(1*(q[(n*(t0-1)+1):(n*t)]>gamma_hat)-1*(q[(n*(t0-2)+1):(n*(t-1))]>gamma_hat),
              x2[(n*(t0-1)+1):(n*t),]*(q[(n*(t0-1)+1):(n*t)]>gamma_hat)-x2[(n*(t0-2)+1):(n*(t-1)),]*(q[(n*(t0-2)+1):(n*(t-1))]>gamma_hat))
  
  
  X <- cbind(X1,X2)
  resid <- Y-X%*%alpha_hat
  
  g_bar <- t(Z) %*% resid/n
  
  # 2. Resampling
  
  if(is.null(rng)){
    i_star <- sample(n,n,replace=TRUE)
  }else{
    i_star <- rng
  }
  i_star_T_large <- rep(NA, length=n*t)
  i_star_T_small <- rep(NA, length=n*(t-2))
  for(i in 1:t){
    i_star_T_large[((i-1)*n+1):(i*n)] <- c(i_star + (i-1)*n)
  }
  for(i in 1:(t-(t0-1))){
    i_star_T_small[((i-1)*n+1):(i*n)] <- c(i_star + (i-1)*n)
  }
  
  # 3. Generate DATA
  # - y_star
  # - Y_star (by residual/null imposed)
  # - X_star
  # - new_data (needed for constructing regressors)
  
  Z_star <- Z[i_star_T_small,]
  y_star <- y[i_star_T_large]
  q_star <- q[i_star_T_large]
  x1_star <- x1[i_star_T_large,]
  x2_star <- x2[i_star_T_large,]
  
  X1 <- FDx[(n*(t0-1)+1):(n*t),]
  X2 <- cbind(1*(q[(n*(t0-1)+1):(n*t)]>gamma)-1*(q[(n*(t0-2)+1):(n*(t-1))]>gamma),
              x2[(n*(t0-1)+1):(n*t),]*(q[(n*(t0-1)+1):(n*t)]>gamma)-x2[(n*(t0-2)+1):(n*(t-1)),]*(q[(n*(t0-2)+1):(n*(t-1))]>gamma))
  
  X <- cbind(X1,X2)
  
  Y_star <- X[i_star_T_small,]%*%alpha+resid[i_star_T_small]
  
  return(list(Z=Z_star, Y=Y_star, y=y_star, x1=x1_star, x2=x2_star, q=q_star, g_bar=g_bar))
}

#=================================================#
#     Unconstrained estimation for bootstrap      #
#=================================================#

# For fixed gamma
boot_fixed_est <- function(Y, 
                           x1, 
                           x2, 
                           q, 
                           gamma, 
                           g_bar0, 
                           weight=NULL, 
                           weight_out=FALSE, 
                           inst, 
                           n=NULL, 
                           t=NULL,
                           t0=NULL){
  
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
  FDx <- matrix(NA, nrow=n*t, ncol=p1)
  for(i in 1:p1){
    FDx[,i] <-FD_transform(x1[,i], n, t)
  }
  
  X1 <- FDx[(n*(t0-1)+1):(n*t),]
  X2 <- cbind(1*(q[(n*(t0-1)+1):(n*t)]>gamma)-1*(q[(n*(t0-2)+1):(n*(t-1))]>gamma),
              x2[(n*(t0-1)+1):(n*t),]*(q[(n*(t0-1)+1):(n*t)]>gamma)-x2[(n*(t0-2)+1):(n*(t-1)),]*(q[(n*(t0-2)+1):(n*(t-1))]>gamma))
  
  X <- cbind(X1,X2)
  Z <- inst
  
  if(is.null(weight)){
    W <- diag(rep(1,ncol(Z)))
  }else{
    W <- weight
  }
  
  g1_bar <- t(Z) %*% (Y)/n
  g2_bar <- t(Z) %*% (X)/n
  
  singular <- FALSE
  tryCatch(  
    coef <- solve(t(g2_bar) %*% W %*% g2_bar, t(g2_bar) %*% W %*% (g1_bar-g_bar0)),
    error=function(e){singular <<- TRUE}
  )
  
  if(singular==FALSE){
    g_bar <- g1_bar - g2_bar %*% coef - g_bar0   # GMM objective
    obj <- t(g_bar) %*% W %*% g_bar     # GMM objective
    
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


# For gamma grid
boot_unconstrained_est <- function(Y,
                                   x1, 
                                   x2, 
                                   q, 
                                   q_grid, 
                                   g_bar0,
                                   two_step=TRUE, 
                                   weight=NULL, 
                                   weight_out=FALSE, 
                                   instmethod=1, 
                                   inst=NULL,
                                   n=NULL, 
                                   t=NULL,
                                   t0=3){
  
  if(instmethod==1){
    inst <- inst
    objlist <- rep(NA, length(q_grid))
    
    # First step
    for(q_iter in 1:length(q_grid)){
      boot0 <- boot_fixed_est(Y, x1, x2, q, q_grid[q_iter], g_bar0, weight, weight_out=F, inst, n, t, t0)
      objlist[q_iter] <- boot0$objective
    }
    
    # Second Step
    if(two_step==TRUE){
      ix <- which.min(objlist)
      weight <- boot_fixed_est(Y, x1, x2, q, q_grid[ix], g_bar0, weight, weight_out=T, inst, n, t, t0)$weight
      for(q_iter in 1:length(q_grid)){
        boot0 <- boot_fixed_est(Y, x1, x2, q, q_grid[q_iter], g_bar0, weight, weight_out=F, inst, n, t, t0)
        objlist[q_iter] <- boot0$objective
      }
    }
    
    ix <- which.min(objlist)
    boot <- boot_fixed_est(Y, x1, x2, q, q_grid[ix], g_bar0, weight, weight_out=weight_out, inst, n, t, t0)
    if(boot$objective!=objlist[ix]) stop()
  }
  
  return(list(coefficients=boot$coefficients, 
              threshold=q_grid[ix], 
              objective=boot$objective, 
              objlist=objlist, 
              used_grid=q_grid, 
              weight=boot$weight,
              g_bar=boot$g_bar))
}





#=================================================#
#      Constrained estimation for bootstrap       #
#=================================================#
# - Continuity imposed


# for fixed gamma
boot_fixed_constrained_est <- function(Y, 
                                       x1, 
                                       q, 
                                       gamma, 
                                       g_bar0, 
                                       weight=NULL, 
                                       weight_out=FALSE, 
                                       inst, 
                                       n=NULL, 
                                       t=NULL,
                                       t0=NULL){
  
  if((is.null(n)|is.null(t))){
    stop("The dimension of panel data should be determined either by providing 'n' and 't'")
  }
  if(is.null(inst)){
    stop("'inst' is not given.")
  }
  
  p1 <- length(x1)/(n*t)
  
  x1 <- as.matrix(x1)
  
  FDx <- matrix(NA, nrow=n*t, ncol=p1)
  for(i in 1:p1){
    FDx[,i] <-FD_transform(x1[,i], n, t)
  }
  
  X1 <- FDx[(n*(t0-1)+1):(n*t),]
  X2 <- (q[(n*(t0-1)+1):(n*t)]-gamma)*(q[(n*(t0-1)+1):(n*t)]>gamma)-(q[(n*(t0-2)+1):(n*(t-1))]-gamma)*(q[(n*(t0-2)+1):(n*(t-1))]>gamma)
  
  X <- cbind(X1,X2)
  Z <- inst
  
  if(is.null(weight)){
    W <- diag(rep(1,ncol(Z)))
  }else{
    W <- weight
  }
  
  g1_bar <- t(Z) %*% (Y)/n
  g2_bar <- t(Z) %*% (X)/n
  
  singular <- FALSE
  tryCatch(  
    coef <- solve(t(g2_bar) %*% W %*% g2_bar, t(g2_bar) %*% W %*% (g1_bar-g_bar0)),
    error=function(e){singular <<- TRUE}
  )
  
  if(singular==FALSE){
    g_bar <- g1_bar - g2_bar %*% coef - g_bar0   # GMM objective
    obj <- t(g_bar) %*% W %*% g_bar     # GMM objective
    
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

# for gamma grid
boot_constrained_est <- function(Y,
                                 x1, 
                                 q, 
                                 q_grid, 
                                 g_bar0,
                                 two_step=TRUE, 
                                 weight=NULL, 
                                 weight_out=FALSE, 
                                 instmethod=1, 
                                 inst=NULL,
                                 n=NULL, 
                                 t=NULL,
                                 t0=3){
  
  
  if(instmethod==1){
    inst <- inst
    
    objlist <- rep(NA, length(q_grid))
    
    # First step
    for(q_iter in 1:length(q_grid)){
      boot0 <- boot_fixed_constrained_est(Y, x1, q, q_grid[q_iter], g_bar0, weight, weight_out=F, inst, n, t, t0)
      objlist[q_iter] <- boot0$objective
    }
    
    # Second Step
    if(two_step==TRUE){
      ix <- which.min(objlist)
      weight <- boot_fixed_constrained_est(Y, x1, q, q_grid[ix], g_bar0, weight, weight_out=T, inst, n, t, t0)$weight
      for(q_iter in 1:length(q_grid)){
        boot0 <- boot_fixed_constrained_est(Y, x1, q, q_grid[q_iter], g_bar0, weight, weight_out=F, inst, n, t, t0)
        objlist[q_iter] <- boot0$objective
      }
    }
    
    ix <- which.min(objlist)
    boot <- boot_fixed_constrained_est(Y, x1, q, q_grid[ix], g_bar0, weight, weight_out=weight_out, inst, n, t, t0)
  }
  
  return(list(coefficients=boot$coefficients, 
              threshold=q_grid[ix], 
              objective=boot$objective, 
              objlist=objlist, 
              used_grid=q_grid, 
              weight=boot$weight,
              g_bar=boot$g_bar))
}