
rm(list=ls())

library(readxl)
library(doRNG)
library(parallel)
library(foreach)
library(doParallel)

source("ftn_est.R")
source("ftn_boot.R")


nc <- detectCores()  # Number of clusters

filenames <- c("LEV_lin.Rdata", "TQ_lin.Rdata", "TQ_submin_lin.Rdata")

for(iter in c(1,2,3)){
  
  # setting parameter
  
  n = 1459
  t = 10
  B = 1000
  
  raw <- read_excel("data1.xls")
  dat <- matrix(nrow=n*t, ncol=5)
  tq <- lev <- matrix(nrow=n*t, ncol=1)
  pdat <- ncol(dat)
  for(i in 1:t){
    if(iter==1){
      dat[(n*(i-1)+1):(n*i),] <- cbind(as.matrix(raw[seq(1, n*(t+1), by=(t+1))+i,c("inv", "cf")]),
                                       as.matrix(raw[seq(1, n*(t+1), by=(t+1))+(i-1),c("ppe", "roa", "lev")]))
    }else{
      dat[(n*(i-1)+1):(n*i),] <- cbind(as.matrix(raw[seq(1, n*(t+1), by=(t+1))+i,c("inv","cf")]),
                                       as.matrix(raw[seq(1, n*(t+1), by=(t+1))+(i-1),c("ppe", "roa", "tq")]))
    }
    lev[(n*(i-1)+1):(n*i),1] <- as.matrix(raw[seq(1, n*(t+1), by=(t+1))+(i-1),c("lev")])
    tq[(n*(i-1)+1):(n*i),1] <- as.matrix(raw[seq(1, n*(t+1), by=(t+1))+(i-1),c("tq")])
  }
  if(iter==1){
    names <- c("lagI", "CF", "PPE", "ROA", "LEV")
    names2 <- c("lagI", "CF", "PPE", "ROA", "LEV")
    
    y <- dat[,1]
    x <- cbind(lag_dpd(y,n,t,1), dat[,-c(1)])
    xx <- x[,which(names %in% names2)]
  }else if(iter==2){
    names <- c("lagI", "CF", "PPE", "ROA", "TQ")
    names2 <- c("lagI", "CF", "PPE", "ROA", "TQ")
    
    y <- dat[,1]
    x <- cbind(lag_dpd(y,n,t,1), dat[,-c(1)])
    xx <- as.matrix(x[,which(names %in% names2)])
  }else{
    names <- c("lagI", "CF", "PPE", "ROA", "TQ")
    names2 <- c("TQ")
    
    y <- dat[,1]
    x <- cbind(lag_dpd(y,n,t,1), dat[,-c(1)])
    xx <- as.matrix(x[,which(names %in% names2)])
  }
  
  ix <- rep(TRUE,n)
  nix <- rep(0,n)
  
  for(i in 1:t){
    for(j in 1:(pdat-1)){
      ix <- ix & ifelse(is.na(!(dat[(n*(i-1)+1):(n*i),j]<quantile(dat[,j],0.005,na.rm=TRUE)|
                                  dat[(n*(i-1)+1):(n*i),j]>quantile(dat[,j],0.995,na.rm=TRUE))),
                        TRUE,
                        !(dat[(n*(i-1)+1):(n*i),j]<quantile(dat[,j],0.005,na.rm=TRUE)|
                            dat[(n*(i-1)+1):(n*i),j]>quantile(dat[,j],0.995,na.rm=TRUE)))
    }
    if(iter==1){
      nix <- nix + ifelse(lev[(n*(i-1)+1):(n*i),1] > 0, 1, 0)
    }
    if(iter!=1){
      nix <- nix + ifelse(tq[(n*(i-1)+1):(n*i),1] < 5 , 1, 0)
    }
  }
  ix <- ix & (nix>5)
    
  n <- sum(ix)
  y <- y[which(rep(ix,t)==TRUE)]
  x <- x[which(rep(ix,t)==TRUE),]
  xx <- as.matrix(xx[which(rep(ix,t)==TRUE),])
  colnames(x) <- names
  colnames(xx) <- names2
  q <- x[,pdat]
  rm(raw,dat)
  
  p1 <- ncol(x)
  p2 <- ncol(xx)
  p <- p1+p2+1
  
  inst <- arrange_inst(cbind(lag_dpd_inst(y,n,t,lb=2,ub=2,t0=3),
                             lag_dpd_inst(x[,2],n,t,lb=1,ub=1,t0=3),
                             lag_dpd_inst(x[,3],n,t,lb=1,ub=1,t0=3),
                             lag_dpd_inst(x[,4],n,t,lb=1,ub=1,t0=3),
                             lag_dpd_inst(x[,5],n,t,lb=1,ub=1,t0=3)),
                       n, t, t0=3)
  
  q_grid <- quantile(q, seq(0.1,0.9,0.01))
  
  gmm <- gmm_unconstrained_est(y, x, xx, q, q_grid, n=n, t=t, t0=3, inst=inst) 
  alpha_hat <- gmm$coefficients 
  gamma_hat <- gmm$threshold
  
  # Linearity Test Statistics
  Wald <- rep(NA, length(q_grid))
  
  for(i in 1:length(q_grid)){
    W_i <- gmm_unconstrained_est(y, x, xx, q, q_grid[i], inst=inst, n=n, t=t, weight_out=TRUE, two_step=FALSE)$weight
    gmm_i <- gmm_unconstrained_est(y, x, xx, q, q_grid[i], inst=inst, n=n, t=t, weight_out=TRUE, two_step=FALSE)
    delta_i <- gmm_i$coefficients[-c(1:p1)]
    Omega <- solve(gmm_i$weight)
    Mn_i <- Design_est(inst, y, x, xx, q, q_grid[i], delta_i, delta_i[p2+1], n=n, t=t)$Mn
    Wald[i] <- n*t(delta_i) %*% solve((solve(t(Mn_i)%*%W_i%*%Mn_i)%*%t(Mn_i)%*%W_i%*%Omega%*%W_i%*%Mn_i%*%solve(t(Mn_i)%*%W_i%*%Mn_i))[(p1+1):(p1+p2+1),(p1+1):(p1+p2+1)]) %*% delta_i
  }
  supW <- max(Wald)
  
  # Linearity Test Bootstrap
  alpha_null <- rep(0, p)
  alpha_null[1:p1] <- alpha_hat[1:p1]
  
  cl <- makeCluster(nc)
  registerDoParallel(cl)
  registerDoRNG(200) #RNG

  try({
    supW_boot <- foreach(b=1:B, .combine=c, .errorhandling="pass") %dopar% {
      istar <- sample(n,n,replace=TRUE)
      
      booted <- boot_null(inst, y, x, xx, q, alpha_null, alpha_hat, 0, gamma_hat, n=n, t=t, t0=3, rng=istar)
      Z_star <- booted$Z
      Y_star <- booted$Y
      y_star <- booted$y
      x1_star <- booted$x1
      x2_star <- booted$x2
      q_star <- booted$q
      g_bar <- booted$g_bar
      
      Wald_i <- rep(NA, length(q_grid))
      
      for(i in 1:length(q_grid)){
        W_istar <- boot_unconstrained_est(Y_star, x1_star, x2_star, q_star, q_grid[i], g_bar=g_bar, inst=Z_star, n=n, t=t, t0=3, weight_out=TRUE, two_step=FALSE)$weight
        boot_i <- boot_unconstrained_est(Y_star, x1_star, x2_star, q_star, q_grid[i], g_bar=g_bar, inst=Z_star, n=n, t=t, t0=3, weight=W_istar, weight_out=TRUE, two_step=FALSE)
        delta_istar <- boot_i$coefficients[-c(1:p1)]
        Omega_star <- solve(boot_i$weight)
        Mn_istar <- Design_est(Z_star, y_star, x1_star, x2_star, q_star, q_grid[i], delta_istar, delta_istar[p2+1], n=n, t=t)$Mn
        Wald_i[i] <- n*t(delta_istar)%*% solve((solve(t(Mn_istar)%*%W_istar%*%Mn_istar)%*%t(Mn_istar)%*%W_istar%*%Omega_star%*%W_istar%*%Mn_istar%*%solve(t(Mn_istar)%*%W_istar%*%Mn_istar))[(p1+1):(p1+p2+1),(p1+1):(p1+p2+1)]) %*% delta_istar
      }
      
      return(max(Wald_i))
    }
    
    LT <- mean(1*(c(supW) < supW_boot))
    print(LT)
  })
  stopCluster(cl)
  
  save.image(filenames[iter])
}