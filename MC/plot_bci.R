library(np)

load("./boot_n1600_g0.25_p0.Rdata")



# GAMMA



x0 <- theta_hat[,3]-delta[1]
xu <- cvs_naive[,3]
ix <- (abs(x0)<quantile(abs(x0),.9))

pdf("plot_bci1.pdf")
par(mar=c(5,5.5,1,1))
plot(x0[ix], xu[ix], ylim=c(0,4), xlab=expression(hat(delta)[1]-delta[10]), ylab=expression(paste(".95th quantile of  ", abs(hat(delta)[1]^{"*"}-hat(delta)[1]))),
     cex.lab=1.5, pch=15)
lines(seq(-2,2,0.01),abs(seq(-2,2,0.01)), col="red")
dev.off()


xu <- cvs_perc2[,3]

pdf("plot_bci2.pdf")
par(mar=c(5,5.5,1,1))
plot(x0[ix], xu[ix], ylim=c(0,4), xlab=expression(hat(delta)[1]-delta[10]), ylab=expression(paste(".95th quantile of  ", abs(hat(delta)[1]^{"*"}-delta[10]^{"*"}))), 
     cex.lab=1.5, pch=15)
lines(seq(-2,2,0.01),abs(seq(-2,2,0.01)), col="red")
dev.off()




# DELTA 3


x0 <- theta_hat[,3]-delta[1]
xu <- cvs_naive_l[,3]
xl <- cvs_naive_u[,3]
ix <- (x0<quantile(x0,.95))&(x0>quantile(x0,.05))


pdf("plot_bci3.pdf")
par(mar=c(5,5.5,1,1))
plot(x0[ix], xu[ix], ylim=c(-4,4), xlab=expression(hat(delta)[1]-delta[10]), ylab=expression(paste(".025/.975th quantiles of  ", (hat(delta)[1]^{"*"}-hat(delta)[1]))),
     cex.lab=1.5, pch=15, col="blue")
par(new=TRUE)
plot(x0[ix], xl[ix], ylim=c(-4,4), ylab="", xlab="", pch=16)
abline(h=0,col='black',lty=2)
abline(0,1,col='red')
legend(x = "topright", legend=c("97.5th boot-quantile", "2.5th boot-quantile"), pch=c(15,16), col = c("blue","black")) 
dev.off()



xu <- cvs_perc2_l[,3]
xl <- cvs_perc2_u[,3]


pdf("plot_bci4.pdf")
par(mar=c(5,5.5,1,1))
plot(x0[ix], xu[ix], ylim=c(-4,4), xlab=expression(hat(delta)[1]-delta[10]), ylab=expression(paste(".025/.975th quantile of  ", (hat(delta)[1]^{"*"}-delta[10]^{"*"}))), 
     cex.lab=1.5, pch=15, col="blue")
par(new=TRUE)
plot(x0[ix], xl[ix], ylim=c(-4,4), ylab="", xlab="", pch=16)
abline(h=0,col='black',lty=2)
abline(0,1,col='red')
legend(x = "topright", legend=c("97.5th boot-quantile", "2.5th boot-quantile"), pch=c(15,16), col = c("blue","black")) 
dev.off()