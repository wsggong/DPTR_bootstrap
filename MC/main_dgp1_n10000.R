rm(list=ls())

source("ftn_dgp1.R")
source("ftn_est.R")
source("ftn_boot.R")

library(doRNG)
library(parallel)
library(foreach)
library(doParallel)

nc <- detectCores()/2   # Parallelization number of cores

# setting parameter

M <- 1000                          # Number of MC repetitions
B1 <- 500
B <- 500                           # Number of bootstrap iterations
n_list <- c(10000)  # Sample size
gamma_grid <-  c(0.25)             # Threshold location
delta1_grid <- c(0, 1)        # Jump size
t <- 6

beta <- c(0, 0.6, 1)
delta <- c(NA, 0, 2)
sigma <- 0.5
p1 <- 2
p2 <- 2
p <- p1+p2+1
pp <- p+1

for(param_iter1 in 1:length(gamma_grid)){
  for(param_iter2 in 1:length(delta1_grid)){
    for(n_iter in 1:length(n_list)){
    cat("\n")  
    n <- n_list[n_iter]
    gamma <- gamma_grid[param_iter1]
    delta[1] <- delta1_grid[param_iter2]-delta[p2+1]*gamma
    cat("n=",n, "gamma=", gamma, "delta1=", delta[1], "\n")
    
    param0 <- c(beta[-1],delta,gamma)
    
    t1 <- Sys.time()
    
    cl <- makeCluster(nc)
    registerDoParallel(cl)
    registerDoRNG(200) #RNG
    
    res <- foreach(i=1:M, .combine=rbind, .errorhandling="pass") %dopar% {
      dat <- dgp(n, t, beta, delta, gamma, sigma)
      y <- dat[,3]
      q <- dat[,4]
      x <- cbind(lag_dpd(y,n,t,1),q)
      q_grid <- sort(c(gamma,quantile(dat[,4], seq(.05,.95,.02))))
      inst <- arrange_inst(cbind(lag_dpd_inst(y,n,t,lb=2,t0=3),
                                 lag_dpd_inst(q,n,t,lb=1,t0=3)),
                           n, t, t0=3)
      
      gmm <- gmm_unconstrained_est(y, x, x, q, q_grid, n=n, t=t, t0=3, inst=inst) 
      alpha_hat <- gmm$coefficients 
      gamma_hat <- gmm$threshold
      theta_hat <- c(alpha_hat, gamma_hat)
      Wn <- gmm$weight
      Ln_hat <- gmm$objective
      
      # Continuity Test Statistics
      
      con_gmm <- gmm_constrained_est(y, x, q, q_grid, n=n, t=t, t0=3, inst=inst, weight=Wn, two_step=FALSE)
      Ln_tilde <- con_gmm$objective
      LR <- n*(Ln_tilde-Ln_hat)

      # Continuity Test Bootstrap

      gamma_kink <- con_gmm$threshold
      alpha_kink <- rep(0, p)
      alpha_kink[1:p1] <- con_gmm$coefficients[1:p1]
      alpha_kink[p1+1] <- -gamma_kink*con_gmm$coefficients[p1+1]
      alpha_kink[p] <- con_gmm$coefficients[p1+1]
      
      # Continuity Test Bootstrap
      
      gamma_kink <- con_gmm$threshold
      alpha_kink <- rep(0, p)
      alpha_kink[1:p1] <- con_gmm$coefficients[1:p1]
      alpha_kink[p1+1] <- -gamma_kink*con_gmm$coefficients[p1+1]
      alpha_kink[p] <- con_gmm$coefficients[p1+1]
      
      LRsims <- rep(NA, B1)
      for(b in 1:B1){
        istar <- sample(n,n,replace=TRUE)
        
        booted <- boot_null(inst, y, x, x, q, alpha_kink, alpha_hat, gamma_kink, gamma_hat, n=n, t=t, t0=3, rng=istar)
        Z_star <- booted$Z
        Y_star <- booted$Y
        x1_star <- booted$x1
        x2_star <- booted$x2
        q_star <- booted$q
        g_bar <- booted$g_bar
        
        uncon_boot_est <- boot_unconstrained_est(Y_star, x1_star, x2_star, q_star, q_grid, g_bar0=g_bar, n=n, t=t, t0=3, inst=Z_star)
        Ln_hat_star <- uncon_boot_est$objective
        Wn_star <- uncon_boot_est$weight
        con_boot_est <- boot_constrained_est(Y_star, x1_star, q_star, q_grid, g_bar0=g_bar, two_step=FALSE, weight=Wn_star, n=n, t=t, t0=3, inst=Z_star)
        Ln_tilde_star <- con_boot_est$objective
        
        LRsims[b] <- n*(Ln_tilde_star-Ln_hat_star)
      }
      
      Chat2 <- quantile(LRsims, .50)
      
      # BOOTSTRAP
      
      boot_res <- matrix(0, nrow=B, ncol=(2*pp))
      
      # BOOTSTRAP Population parameter
      Ratio <- min(LR/(Chat2*n^(1/4)),1)
      alpha_null1 <- Ratio*alpha_hat+(1-Ratio)*alpha_kink
      gamma_null1 <- Ratio*gamma_hat+(1-Ratio)*gamma_kink
      
      
      for(b in 1:B){
        istar <- sample(n,n,replace=TRUE)
        
        # R-B
        
        booted <- boot_null(inst, y, x, x, q, alpha_null1, alpha_hat, gamma_null1, gamma_hat, n=n, t=t, rng=istar)
        Z_star <- booted$Z
        Y_star <- booted$Y
        x1_star <- booted$x1
        x2_star <- booted$x2
        q_star <- booted$q
        g_bar <- booted$g_bar
        
        uncon_boot_est <- boot_unconstrained_est(Y_star, x1_star, x2_star, q_star, q_grid-gamma+gamma_null1, g_bar0=g_bar, n=n, t=t, inst=Z_star)
        
        boot_res[b,(1):(pp)] <- c(uncon_boot_est$coefficients, uncon_boot_est$threshold) - c(alpha_null1, gamma_null1)
        
        
        # Naive Bootstrap
        
        booted <- boot_null(inst, y, x, x, q, alpha_hat, alpha_hat, gamma_hat, gamma_hat, n=n, t=t, rng=istar)
        Z_star <- booted$Z
        Y_star <- booted$Y
        x1_star <- booted$x1
        x2_star <- booted$x2
        q_star <- booted$q
        g_bar <- booted$g_bar
        
        uncon_boot_est <- boot_unconstrained_est(Y_star, x1_star, x2_star, q_star, q_grid-gamma+gamma_hat, g_bar0=g_bar, n=n, t=t, inst=Z_star)
        
        boot_res[b,(pp+1):(2*pp)] <- c(uncon_boot_est$coefficients, uncon_boot_est$threshold) - c(alpha_hat, gamma_hat)
        
        
      }
      
      # Asymptotic CI
      Des <- Design_est(inst, y, x, x, q, gamma_hat, delta=alpha_hat[(p1+1):p], delta3=alpha_hat[p], n=n, t=t, t0=3)
      Mn <- Des$Mn
      Gn <- Des$Gn
      Hn <- Des$Hn
      
      Omega_inv <- gmm_unconstrained_est(y, x, x, q, gamma_hat, weight=Wn, weight_out=TRUE, two_step=FALSE, n=n, t=t, inst=inst)$weight
      
      A <- cbind(Mn, Gn)
      sd_hat <- sqrt(diag(solve(t(A) %*% Omega_inv %*% A)/n))
      
      Ve_inv <- t(Mn)%*%Omega_inv%*%Mn
      Le <- eigen(Ve_inv)$values
      Pe <- eigen(Ve_inv)$vectors
      Ve_sqrt <- Pe%*%diag(1/sqrt(Le))%*%t(Pe)
      Ve <- Pe%*%diag(1/Le)%*%t(Pe)
      
      e <- Ve_sqrt %*% matrix(rnorm(1000*ncol(Mn)), nrow=ncol(Mn), ncol=1000)
      
      Vv_meat <- Omega_inv-Omega_inv%*%(Mn)%*%Ve%*%t(Mn)%*%Omega_inv
      Vv_inv <- t(Hn)%*%Vv_meat%*%Hn
      Vv_sqrt <- (1/sqrt(Vv_inv))
      
      v <- matrix(pmax(c(Vv_sqrt)*rnorm(1000), 0), nrow=1, ncol=1000)
      
      sim <- e - Ve%*%t(Mn)%*%Omega_inv%*%Hn%*%v
      
      cvs_d <- sd_hat*qnorm(.975)
      cvs_c <- c(apply(abs(sim), 1, quantile, .95)/sqrt(n), sqrt(quantile(v, .95)/sqrt(n)))
      
      return(c(theta_hat, 
               apply(abs(boot_res[,1:(2*pp)]),2,quantile,.95),
               apply((boot_res[,1:(2*pp)]),2,quantile,.975),
               apply((boot_res[,1:(2*pp)]),2,quantile,.025),
               cvs_c, cvs_d, LR))
    }
    stopCluster(cl)
    
    t2 <- Sys.time()
    print(t2-t1)
    
    
    theta_hat <- res[,1:(pp)]
    cvs_perc2 <- res[,(pp+1):(2*pp)]
    cvs_naive <- res[,(2*pp+1):(3*pp)]
    
    cvs_perc2_l <- res[,(3*pp+1):(4*pp)]
    cvs_naive_l <- res[,(4*pp+1):(5*pp)]
    
    cvs_perc2_u <- res[,(5*pp+1):(6*pp)]
    cvs_naive_u <- res[,(6*pp+1):(7*pp)]
    
    acvs_c <- res[,(7*pp+1):(8*pp)]
    acvs_d <- res[,(8*pp+1):(9*pp)]
    
    param0_mat <- matrix(rep(param0,M), nrow=M, ncol=p+1, byrow=T)
    test_param <- theta_hat-param0_mat
    
    cr95_2 <- colMeans(abs(test_param)<cvs_perc2)
    cr95_4 <- colMeans(abs(test_param)<cvs_naive)
    
    cr95a_2 <- colMeans(test_param<cvs_perc2_l & test_param>cvs_perc2_u)
    cr95a_4 <- colMeans(test_param<cvs_naive_l & test_param>cvs_naive_u)
    
    
    acr95_1 <- colMeans(abs(test_param)<acvs_c)
    acr95_2 <- colMeans(abs(test_param)<acvs_d)
    
    bci2_l <- round(apply(cvs_perc2,2,median),3)
    npbci_l <- round(apply(cvs_naive,2,median),3)
    aci_l <- round(apply(acvs_d,2,median),3)
    
    filename <- paste("boot_n",n,"_g",gamma,"_p",delta1_grid[param_iter2],"2.Rdata",sep="")
    save.image(filename)
    
    cat("===========================================================\n", file="text_output.txt", append=TRUE)
    cat(paste(as.character(Sys.time()),"\n"), file="text_output.txt", append=TRUE)
    cat(paste("filename: ", filename, "\n", sep="  "), file="text_output.txt", append=TRUE)
    cat(paste("sd  : ", round(sd(theta_hat[,1]),3), round(sd(theta_hat[,2]),3), round(sd(theta_hat[,3]),3), round(sd(theta_hat[,4]),3),
              round(sd(theta_hat[,5]),3),round(sd(theta_hat[,6]),3), "\n", sep="  "), file="text_output.txt", append=TRUE)
    cat(paste("RB2","\n", sep="  "), file="text_output.txt", append=TRUE)
    cat(paste("cr95: ", cr95a_2[1], cr95a_2[2], cr95a_2[3], cr95a_2[4], cr95a_2[5], cr95a_2[6], "\n", sep="  "), file="text_output.txt", append=TRUE)
    cat(paste("NP-B","\n", sep="  "), file="text_output.txt", append=TRUE)
    cat(paste("cr95: ", cr95a_4[1], cr95a_4[2], cr95a_4[3], cr95a_4[4], cr95a_4[5], cr95a_4[6], "\n", sep="  "), file="text_output.txt", append=TRUE)
    cat(paste("Asymptotic CI (Continuous)","\n", sep="  "), file="text_output.txt", append=TRUE)
    cat(paste("cr95: ", acr95_1[1], acr95_1[2], acr95_1[3], acr95_1[4], acr95_1[5], acr95_1[6], "\n", sep="  "), file="text_output.txt", append=TRUE)
    cat(paste("Asymptotic CI (Discontinuous)","\n", sep="  "), file="text_output.txt", append=TRUE)
    cat(paste("cr95: ", acr95_2[1], acr95_2[2], acr95_2[3], acr95_2[4], acr95_2[5], acr95_2[6], "\n", sep="  "), file="text_output.txt", append=TRUE)
    cat(paste("Elapsed time: ", round(difftime(t2,t1, units=c("mins")),5), "mins", "\n"), file="text_output.txt", append=TRUE)
    }
  }
}
