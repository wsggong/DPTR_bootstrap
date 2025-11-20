dgp <- function(n,t,beta,delta,gamma,sigma){ # weakly exogenous dgp
  
  preheating <- 20
  y <- matrix(NA, nrow=n, ncol=t+preheating)
  x1 <- matrix(NA, nrow=n, ncol=t+preheating)
  d <- matrix(NA, nrow=n, ncol=t+preheating)
  res <- matrix(NA, nrow=n*t, ncol=4)
  
  y[,1] <- rnorm(n)
  x1[,1] <- rnorm(n)
  oe <- rnorm(n)
  for(i in 2:(t+preheating)){
    e <- rnorm(n)
    v <- rnorm(n)
    x1[,i] <- 0.7*x1[,i-1]+(v+oe)/sqrt(2)
    y[,i] <- beta[1]+beta[2]*y[,i-1]+beta[3]*x1[,i]+
      ifelse(x1[,i] <= gamma, 0, delta[1]+delta[2]*y[,i-1]+delta[3]*x1[,i])+sigma*(e)
    oe <- e
  }
  
  res <- cbind(rep(1:n,t), 
               rep(1:t,each=n),
               c(y[,c((preheating+1):(preheating+t))]),
               c(x1[,c((preheating+1):(preheating+t))]))
  return(res)
}