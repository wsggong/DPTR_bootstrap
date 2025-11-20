rm(list=ls())

source("ftn_dgp1.R")
source("ftn_est.R")
source("ftn_boot.R")

library(doRNG)
library(parallel)
library(foreach)
library(doParallel)

nc <- detectCores()   # Parallelization number of cores

# setting parameter

M <- 2000                           # Number of MC repetitions
B <- 500                           # Number of bootstrap iterations
n_list <- c(400, 800, 1600)         # Sample size
gamma_grid <-  c(0.25)             # Threshold location
delta1_grid <- c(0, 0.1, 0.2, 0.5, 1)        # Jump size
t <- 6

beta <- c(0, 0.6, 1)
delta <- c(NA, 0, 2)
sigma <- 0.5
p1 <- 2
p2 <- 2
p <- p1+p2+1
pp <- p+1

for(n_iter in 1:length(n_list)){
  for(param_iter1 in 1:length(gamma_grid)){
    for(param_iter2 in 1:length(delta1_grid)){n <- n_list[n_iter];
    cat("\n")
    gamma <- gamma_grid[param_iter1]
    delta[1] <- delta1_grid[param_iter2]-delta[p2+1]*gamma
    cat("n=",n, "gamma=", gamma, "delta1=", delta[1], "\n")
    
    
    param0 <- c(beta[-1],delta,gamma)
    test_grid <- gamma+c(0.1, 0.25, 0.5)
    
    t1 <- Sys.time()
    
    cl <- makeCluster(nc)
    registerDoParallel(cl)
    registerDoRNG(100) #RNG
    
    res <- foreach(i=1:M, .combine=rbind, .errorhandling="pass") %dopar% {
      dat <- dgp(n, t, beta, delta, gamma, sigma)
      y <- dat[,3]
      q <- dat[,4]
      x <- cbind(lag_dpd(y,n,t,1),q)
      
      # Make grid
      q_grid <- sort(c(gamma,quantile(dat[,4], seq(.05,.95,.02))))
      q_grid <- unique(sort(c(test_grid,q_grid)))
      ix <- which(q_grid %in% test_grid)
      
      inst <- arrange_inst(cbind(lag_dpd_inst(y,n,t,lb=2,t0=3),
                                 lag_dpd_inst(dat[,4],n,t,lb=1,t0=3)),
                           n, t, t0=3)
      
      gmm <- gmm_unconstrained_est(y, x, x, q, q_grid, n=n, t=t, t0=3, inst=inst) 
      alpha_hat <- gmm$coefficients 
      gamma_hat <- gmm$threshold
      theta_hat <- c(alpha_hat, gamma_hat)
      Wn <- gmm$weight
      Ln_hat <- gmm$objective
      
      dist <- n*(gmm$objlist-c(gmm$objective))[ix]

      # BOOTSTRAP
      
      boot_res <- matrix(0, nrow=B, ncol=length(test_grid))
      gamma_star <- rep(0, B)
      istar_mat <- matrix(sample(n,B*n,replace=TRUE),nrow=B,ncol=n)
      
      alpha_mat <- matrix(nrow=length(test_grid), ncol=5)
      
      for(test_iter in 1:length(test_grid)){
        gamma_null2 <- test_grid[test_iter]
        con_gmm <- gmm_unconstrained_est(y, x, x, q, gamma_null2, n=n, t=t, inst=inst, weight=Wn, two_step=FALSE)
        alpha_mat[test_iter,] <- con_gmm$coefficients
      }
      
      # Grid Bootstrap
      for(b in 1:B){
        istar <- istar_mat[b,]

        for(test_iter in 1:length(test_grid)){

          alpha_null2 <- alpha_mat[test_iter,]
          gamma_null2 <- test_grid[test_iter]

          booted <- boot_null(inst, y, x, x, q, alpha_null2, alpha_hat, gamma_null2, gamma_hat, n=n, t=t, rng=istar)
          Z_star <- booted$Z
          Y_star <- booted$Y
          x1_star <- booted$x1
          x2_star <- booted$x2
          q_star <- booted$q
          g_bar <- booted$g_bar
          
          uncon_boot_est <- boot_unconstrained_est(Y_star, x1_star, x2_star, q_star, q_grid, g_bar0=g_bar, n=n, t=t, inst=Z_star)
          Ln_hat_star <- uncon_boot_est$objective
          Wn_star <- uncon_boot_est$weight

          con_boot_est <- boot_unconstrained_est(Y_star, x1_star, x2_star, q_star, gamma_null2, g_bar0=g_bar, n=n, t=t, inst=Z_star, weight=Wn_star, two_step=FALSE)
          Ln_tilde_star <- con_boot_est$objective

          boot_res[b,test_iter] <- n*(Ln_tilde_star-Ln_hat_star)

        }
        
        booted <- boot_null(inst, y, x, x, q, alpha_hat, alpha_hat, gamma_hat, gamma_hat, n=n, t=t, rng=istar)
        Z_star <- booted$Z
        Y_star <- booted$Y
        x1_star <- booted$x1
        x2_star <- booted$x2
        q_star <- booted$q
        g_bar <- booted$g_bar
        
        uncon_boot_est <- boot_unconstrained_est(Y_star, x1_star, x2_star, q_star, q_grid, g_bar0=g_bar, n=n, t=t, inst=Z_star)
        gamma_star[b] <- uncon_boot_est$threshold
        
      }
      
      return(c(apply(boot_res,2,quantile,.95), dist, gamma_hat, quantile(abs(gamma_star-gamma_hat),.95)))
    }
    stopCluster(cl)
    
    t2 <- Sys.time()
    print(t2-t1)
    
    cvs_grid <- res[,1:3]
    dist <- res[,4:6]
    gamma_hat <- res[,7]
    cvs_naive <- res[,8]
    
    test_param <- matrix(gamma_hat,nrow=M,ncol=length(test_grid)) - matrix(test_grid,byrow=T,nrow=M,ncol=length(test_grid))
    
    filename <- paste("boot(power)_n",n,"_g",gamma,"_p",delta1_grid[param_iter2],".Rdata",sep="")
    
    power <- round(colMeans(cvs_grid<dist),3)
    power_naive <- round(colMeans(abs(test_param)>matrix(cvs_naive,nrow=M,ncol=length(test_grid))),3)
    
    save.image(filename)
    
    cat("===========================================================\n", file="power_output.txt", append=TRUE)
    cat(paste(as.character(Sys.time()),"\n"), file="power_output.txt", append=TRUE)
    cat(paste("filename: ", filename, "\n", sep="  "), file="power_output.txt", append=TRUE)
    cat(paste("Power","\n", sep="  "), file="power_output.txt", append=TRUE)
    cat(paste("Grid-B:     ", power[1], power[2], power[3], "\n", sep="  "), file="power_output.txt", append=TRUE)
    cat(paste("NP-B(S): ", power_naive[1], power_naive[2], power_naive[3], "\n", sep="  "), file="power_output.txt", append=TRUE)
    
    }
  }
}
