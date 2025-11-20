
rm(list=ls())

library(readxl)
library(doRNG)
library(parallel)
library(foreach)
library(doParallel)

source("ftn_est.R")
source("ftn_boot.R")


nc <- detectCores()  # Number of clusters

thrvar <- c("lev", "tq", "tq")
thrname <- c("LEV", "TQ", "TQ")
filenames <- c("LEV.Rdata", "TQ.Rdata", "TQ_submin.Rdata")

for(iter in 1:3){
  
  # setting constants
  n = 1459
  t = 10
  B1 = 1000 # Bootstrap iterations for continuity test
  B2 = 2000 # Bootstrap iterations for coefficients CIs
  B3 = 500 # Bootstrap iterations for threshold CS
  
  raw <- read_excel("data1.xls")
  dat <- matrix(nrow=n*t, ncol=5)
  
  tq <- lev <- matrix(nrow=n*t, ncol=1)
  pdat <- ncol(dat) # Dimension of regressors
  
  for(i in 1:t){
    if(iter==1){
      # Storing I_t, CF_t, PPE_t-1, ROA_t-1, LEV_t-1 for t=1,...,10
      dat[(n*(i-1)+1):(n*i),] <- cbind(as.matrix(raw[seq(1, n*(t+1), by=(t+1))+i,c("inv", "cf")]),
                                       as.matrix(raw[seq(1, n*(t+1), by=(t+1))+(i-1),c("ppe", "roa", "lev")]))
    }else{
      # Storing I_t, CF_t, PPE_t-1, ROA_t-1, TQ_t-1 for t=1,...,10
      dat[(n*(i-1)+1):(n*i),] <- cbind(as.matrix(raw[seq(1, n*(t+1), by=(t+1))+i,c("inv","cf")]),
                                       as.matrix(raw[seq(1, n*(t+1), by=(t+1))+(i-1),c("ppe", "roa", "tq")]))
    }
    lev[(n*(i-1)+1):(n*i),1] <- as.matrix(raw[seq(1, n*(t+1), by=(t+1))+(i-1),c("lev")])
    tq[(n*(i-1)+1):(n*i),1] <- as.matrix(raw[seq(1, n*(t+1), by=(t+1))+(i-1),c("tq")])
  }
  
  # Defining set of regressors
  if(iter==1){
    # Leverage model (Empirical specification #1 in paper)
    names <- c("lagI", "CF", "PPE", "ROA", "LEV")
    names2 <- c("lagI", "CF", "PPE", "ROA", "LEV")
    
    y <- dat[,1]
    x <- cbind(lag_dpd(y,n,t,1), dat[,-c(1)])
    xx <- x[,which(names %in% names2)]
  }else if(iter==2){
    # TQ model (Empirical specification #2 in paper)
    
    names <- c("lagI", "CF", "PPE", "ROA", "TQ")
    names2 <- c("lagI", "CF", "PPE", "ROA", "TQ")
    
    y <- dat[,1]
    x <- cbind(lag_dpd(y,n,t,1), dat[,-c(1)])
    xx <- as.matrix(x[,which(names %in% names2)])
  }else{
    # restricted TQ model (Empirical specification #3 in paper)
    names <- c("lagI", "CF", "PPE", "ROA", "TQ")
    names2 <- c("TQ")
    
    y <- dat[,1]
    x <- cbind(lag_dpd(y,n,t,1), dat[,-c(1)])
    xx <- as.matrix(x[,which(names %in% names2)])
  }
  
  ix <- rep(TRUE,n) # indexes that will be TRUE if firms remain in sample.
  nix <- rep(0,n) # Storing number of periods that Leverage > 0 or TQ < 5.
  
  for(i in 1:t){
    # For each period, removing observations whose I_t-1, CF_t, PPE_t, ROA_t were not in 0.5%-99.5%
    for(j in 1:(pdat-1)){
      ix <- ix & ifelse(is.na(!(dat[(n*(i-1)+1):(n*i),j]<quantile(dat[,j],0.005,na.rm=TRUE)|
                                  dat[(n*(i-1)+1):(n*i),j]>quantile(dat[,j],0.995,na.rm=TRUE))),
                        TRUE,
                        !(dat[(n*(i-1)+1):(n*i),j]<quantile(dat[,j],0.005,na.rm=TRUE)|
                            dat[(n*(i-1)+1):(n*i),j]>quantile(dat[,j],0.995,na.rm=TRUE)))
    }
    # For leverage threshold model, counting the number of periods that leverage was larger than 0 for each firm
    if(iter==1){
      nix <- nix + ifelse(lev[(n*(i-1)+1):(n*i),1] > 0, 1, 0)
    }
    # For TQ threshold model, counting the number of periods that TQ was smaller than 5 for each firm
    if(iter!=1){
      nix <- nix + ifelse(tq[(n*(i-1)+1):(n*i),1] < 5 , 1, 0)
    }
  }
  
  # For leverage threshold model, we leave firms that had more than five nonzero leverage periods.
  # For TQ threshold model, we leave firms that had more than five periods that TQ was below 5.
  ix <- ix & (nix>5)
  
  n <- sum(ix) # Remaining sample size
  
  y <- y[which(rep(ix,t)==TRUE)] # dependent variable
  x <- x[which(rep(ix,t)==TRUE),] # set of regressors
  xx <- as.matrix(xx[which(rep(ix,t)==TRUE),]) # set of regressors that change coefficients
  q <- x[,pdat] # threshold variable
  
  
  colnames(x) <- names
  colnames(xx) <- names2
  rm(raw,dat)
  
  
  p1 <- ncol(x)
  p2 <- ncol(xx)
  p <- p1+p2+1
  
  # IV for period t is (I_t-2, CF_t-1, PPE_t-2, PPE_t-2, ROA_t-2, TQ/LEV_t-2)
  # t0 = 3, T = 10
  inst <- arrange_inst(cbind(lag_dpd_inst(y,n,t,lb=2,ub=2,t0=3),
                             lag_dpd_inst(x[,2],n,t,lb=1,ub=1,t0=3),
                             lag_dpd_inst(x[,3],n,t,lb=1,ub=1,t0=3),
                             lag_dpd_inst(x[,4],n,t,lb=1,ub=1,t0=3),
                             lag_dpd_inst(x[,5],n,t,lb=1,ub=1,t0=3)),
                       n, t, t0=3)
  
  # Threshold location grid
  q_grid <- quantile(q, seq(0.1,0.9,0.01))
  
  # Unrestricted estimator
  gmm <- gmm_unconstrained_est(y, x, xx, q, q_grid, n=n, t=t, t0=3, inst=inst) 
  alpha_hat <- gmm$coefficients
  gamma_hat <- gmm$threshold
  theta_hat <- c(alpha_hat, gamma_hat)
  Wn <- gmm$weight
  Ln_hat <- gmm$objective
  
  print("Parameter estimator:")
  print(theta_hat)
  plot(q_grid,n*(gmm$objlist-gmm$objective), type="l") # show profiled objective where x-axis is threhsold locations
  
  # Continuity Test Statistics
  con_gmm <- gmm_constrained_est(y, x, q, q_grid, n=n, t=t, t0=3, inst=inst, weight=Wn, two_step=F)
  Ln_tilde <- con_gmm$objective
  LR <- (Ln_tilde-Ln_hat)*n
  
  # Get theta_0^* for continuity Test Bootstrap
  gamma_kink <- con_gmm$threshold
  alpha_kink <- rep(0, p)
  alpha_kink[1:p1] <- con_gmm$coefficients[1:p1]
  alpha_kink[p1+1] <- -gamma_kink*con_gmm$coefficients[p1+1]
  alpha_kink[p] <- con_gmm$coefficients[p1+1]
  
  cl <- makeCluster(nc)
  registerDoParallel(cl)
  registerDoRNG(200) #RNG

  Des <- Design_est(inst, y, x, xx, q, gamma_hat, alpha_hat[(p1+1):p], alpha_hat[p], n, t, t0=3)
  D2 <- cbind(Des$Mn,Des$Hn)
  Matrix::rankMatrix(D2)
  
  try({
    # Bootstrapped continuity test statistics
    LRsims <- foreach(b=1:B1, .combine=c, .errorhandling="pass") %dopar% {
      istar <- sample(n,n,replace=TRUE)

      booted <- boot_null(inst, y, x, xx, q, alpha_kink, alpha_hat, gamma_kink, gamma_hat, n=n, t=t, t0=3, rng=istar)
      Z_star <- booted$Z
      Y_star <- booted$Y
      y_star <- booted$y
      x1_star <- booted$x1
      x2_star <- booted$x2
      q_star <- booted$q
      g_bar <- booted$g_bar

      uncon_boot_est <- boot_unconstrained_est(Y_star, x1_star, x2_star, q_star, q_grid, g_bar0=g_bar, n=n, t=t, t0=3, inst=Z_star)
      Ln_hat_star <- uncon_boot_est$objective
      Wn_star <- uncon_boot_est$weight
      con_boot_est <- boot_constrained_est(Y_star, x1_star, q_star, q_grid, g_bar0=g_bar, two_step=FALSE, weight=Wn_star, n=n, t=t, t0=3, inst=Z_star)
      Ln_tilde_star <- con_boot_est$objective

      return(n*(Ln_tilde_star-Ln_hat_star))
    }


    Chat <- quantile(LRsims,.95) # 5% size critical value for continuity test
    CT <- mean(1*(c(LR) < LRsims)) # p-value of continuity test

    print("p-value for continuity test:")
    print(CT)

    boot_res1 <- matrix(0, nrow=B2, ncol=(p+1))

    # Residual-Bootstrap & Standard Nonparametric Bootstrap

    Ratio <- min(LR/(quantile(LRsims,.5)*n^(1/4)),1)
    alpha_null1 <- Ratio*alpha_hat+(1-Ratio)*alpha_kink
    gamma_null1 <- Ratio*gamma_hat+(1-Ratio)*gamma_kink

    boot_res1 <- foreach(b=1:B2, .combine=rbind, .errorhandling="pass")%dopar%{

      istar <- sample(n,n,replace=TRUE)

      # Residual-Bootstrap

      booted <- boot_null(inst, y, x, xx, q, alpha_null1, alpha_hat, gamma_null1, gamma_hat, n=n, t=t, t0=3, rng=istar)
      Z_star <- booted$Z
      Y_star <- booted$Y
      y_star <- booted$y
      x1_star <- booted$x1
      x2_star <- booted$x2
      q_star <- booted$q
      g_bar <- booted$g_bar

      uncon_boot_est <- boot_unconstrained_est(Y_star, x1_star, x2_star, q_star, q_grid, g_bar0=g_bar, n=n, t=t, t0=3, inst=Z_star)

      return(c(uncon_boot_est$coefficients, uncon_boot_est$threshold) - c(alpha_null1, gamma_null1))
    }

    boot_res11 <- foreach(b=1:B2, .combine=rbind, .errorhandling="pass")%dopar%{

      istar <- sample(n,n,replace=TRUE)

      # Residual-Bootstrap

      booted <- boot_null(inst, y, x, xx, q, alpha_hat, alpha_hat, gamma_hat, gamma_hat, n=n, t=t, t0=3, rng=istar)
      Z_star <- booted$Z
      Y_star <- booted$Y
      y_star <- booted$y
      x1_star <- booted$x1
      x2_star <- booted$x2
      q_star <- booted$q
      g_bar <- booted$g_bar

      uncon_boot_est <- boot_unconstrained_est(Y_star, x1_star, x2_star, q_star, q_grid, g_bar0=g_bar, n=n, t=t, t0=3, inst=Z_star)

      return(c(uncon_boot_est$coefficients, uncon_boot_est$threshold) - c(alpha_hat, gamma_hat))
    }


    cvs_boot <- apply(abs(boot_res1),2,quantile,.95)
    table <- cbind(theta_hat, theta_hat-cvs_boot, theta_hat+cvs_boot)
    rownames(table) <- c(names,paste(c("cons", names2), ".d", sep=""),"tau")

    print("Bootstrap CIs for the coefficients:")
    print(table[-nrow(table),])
  })


  # Grid Bootstrap
  boot_grid <- q_grid

  istar_mat <- matrix(sample(n,n*B3,replace=TRUE), B3, n)

  try({
    boot_res2 <- foreach(qi=1:length(boot_grid), .combine=rbind, .errorhandling="pass")%dopar%{

      gamma <- boot_grid[qi]

      # restricted estimation where the threshold location is gamma
      con_gmm2 <- gmm_unconstrained_est(y, x, xx, q, gamma, n=n, t=t, t0=3, inst=inst, weight=Wn, two_step=FALSE)
      Ln_tilde <- con_gmm2$objective
      alpha_tilde <- con_gmm2$coefficients
      dist2 <- n*(Ln_tilde-Ln_hat) # Distance test statistic

      alpha_null2 <- alpha_tilde
      gamma_null2 <- gamma

      tempboot <- rep(NA,B3)
      for(b in 1:B3){
        istar <- istar_mat[b,]

        booted <- boot_null(inst, y, x, xx, q, alpha_null2, alpha_hat, gamma_null2, gamma_hat, n=n, t=t, t0=3, rng=istar)
        Z_star <- booted$Z
        Y_star <- booted$Y
        y_star <- booted$y
        x1_star <- booted$x1
        x2_star <- booted$x2
        q_star <- booted$q
        g_bar <- booted$g_bar

        uncon_boot_est <- boot_unconstrained_est(Y_star, x1_star, x2_star, q_star, q_grid, g_bar0=g_bar, n=n, t=t, t0=3, inst=Z_star)
        Ln_hat_star <- uncon_boot_est$objective
        Wn_star <- uncon_boot_est$weight
        con_boot_est <- boot_unconstrained_est(Y_star, x1_star, x2_star, q_star, gamma, g_bar0=g_bar, n=n, t=t, t0=3, inst=Z_star, weight=Wn_star, two_step=FALSE)
        Ln_tilde_star <- con_boot_est$objective

        tempboot[b] <- n*(Ln_tilde_star-Ln_hat_star)
      }
      return(c(dist2, quantile(tempboot,c(.9,.95,.99))))
    }
  })
  stopCluster(cl)
  
  save.image(filenames[iter])
}

