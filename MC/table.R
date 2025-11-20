library(xtable)


TAB <- TAB_a <- as.data.frame(matrix(nrow=32+4,ncol=7))
TAB2 <- TAB2_a <- as.data.frame(matrix(nrow=16+4,ncol=6))

TAB[,2] <- c("$I_{t-1}$","",
             "$CF_{t-1}$","",
             "$PPE_{t-1}$","",
             "$ROA_{t-1}$","",
             "$LEV_{t-1}$","",
             "$I_{t-1}$","",
             "$CF_{t-1}$","",
             "$PPE_{t-1}$","",
             "$ROA_{t-1}$","",
             "$LEV_{t-1}$","",
             "intercept","",
             "$I_{t-1}$","",
             "$CF_{t-1}$","",
             "$PPE_{t-1}$","",
             "$ROA_{t-1}$","",
             "$LEV_{t-1}$","",
             "Linearity",
             "Continuity",
             "Threshold",
             "")


TAB[,5] <- c("$I_{t-1}$","",
             "$CF_{t-1}$","",
             "$PPE_{t-1}$","",
             "$ROA_{t-1}$","",
             "$TQ_{t-1}$","",
             "$I_{t-1}$","",
             "$CF_{t-1}$","",
             "$PPE_{t-1}$","",
             "$ROA_{t-1}$","",
             "$TQ_{t-1}$","",
             "intercept","",
             "$I_{t-1}$","",
             "$CF_{t-1}$","",
             "$PPE_{t-1}$","",
             "$ROA_{t-1}$","",
             "$TQ_{t-1}$","",
             "Linearity",
             "Continuity",
             "Threshold",
             "")

TAB2[,2] <- c("$I_{t-1}$","",
             "$CF_{t-1}$","",
             "$PPE_{t-1}$","",
             "$ROA_{t-1}$","",
             "$TQ_{t-1}1\\{TQ_{t-1}\\leq\\gamma\\}$","",
             "$TQ_{t-1}1\\{TQ_{t-1}>\\gamma\\}$","",
             "intercept","",
             "$TQ_{t-1}$","",
             "Linearity",
             "Continuity",
             "Threshold",
             "")

TAB_a[,2] <- c("$I_{t-1}$","",
             "$CF_{t-1}$","",
             "$PPE_{t-1}$","",
             "$ROA_{t-1}$","",
             "$LEV_{t-1}$","",
             "$I_{t-1}$","",
             "$CF_{t-1}$","",
             "$PPE_{t-1}$","",
             "$ROA_{t-1}$","",
             "$LEV_{t-1}$","",
             "intercept","",
             "$I_{t-1}$","",
             "$CF_{t-1}$","",
             "$PPE_{t-1}$","",
             "$ROA_{t-1}$","",
             "$LEV_{t-1}$","",
             "Linearity",
             "Continuity",
             "Threshold",
             "")


TAB_a[,5] <- c("$I_{t-1}$","",
             "$CF_{t-1}$","",
             "$PPE_{t-1}$","",
             "$ROA_{t-1}$","",
             "$TQ_{t-1}$","",
             "$I_{t-1}$","",
             "$CF_{t-1}$","",
             "$PPE_{t-1}$","",
             "$ROA_{t-1}$","",
             "$TQ_{t-1}$","",
             "intercept","",
             "$I_{t-1}$","",
             "$CF_{t-1}$","",
             "$PPE_{t-1}$","",
             "$ROA_{t-1}$","",
             "$TQ_{t-1}$","",
             "Linearity",
             "Continuity",
             "Threshold",
             "")

TAB2_a[,2] <- c("$I_{t-1}$","",
              "$CF_{t-1}$","",
              "$PPE_{t-1}$","",
              "$ROA_{t-1}$","",
              "$TQ_{t-1}1\\{TQ_{t-1}\\leq\\gamma\\}$","",
              "$TQ_{t-1}1\\{TQ_{t-1}>\\gamma\\}$","",
              "intercept","",
              "$TQ_{t-1}$","",
              "Linearity",
              "Continuity",
              "Threshold",
              "")

load("LEV_lin.Rdata")
load("LEV.Rdata")
cvs_boot_l <- apply((boot_res1),2,quantile,0.025)
cvs_boot_u <- apply((boot_res1),2,quantile,0.975)
table <- cbind(theta_hat, theta_hat-cvs_boot, theta_hat+cvs_boot)
table_a <- cbind(theta_hat, theta_hat-cvs_boot_u, theta_hat-cvs_boot_l)
upcvs_boot_l <- apply((boot_res1[,1:p1]+boot_res1[,(p1+2):p]),2,quantile,0.025)
upcvs_boot_u <- apply((boot_res1[,1:p1]+boot_res1[,(p1+2):p]),2,quantile,0.975)
upcvs_boot <- apply(abs(boot_res1[,1:p1]+boot_res1[,(p1+2):p]),2,quantile,0.95)
alpha_up <- alpha_hat[1:p1]+alpha_hat[(p1+2):p]
uptable <- cbind(alpha_up, alpha_up-upcvs_boot, alpha_up+upcvs_boot)
uptable_a <- cbind(alpha_up, alpha_up-upcvs_boot_u, alpha_up-upcvs_boot_l)
TAB[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31),3] <- c(table[1:5,1], uptable[,1], table[6:11,1])
TAB[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32),3] <- c(table[1:5,2], uptable[,2], table[6:11,2])
TAB[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32),4] <- c(table[1:5,3], uptable[,3], table[6:11,3])
TAB[c(33,34),3] <- c(LT, CT)
TAB[35,3] <- gamma_hat
TAB[36,c(3,4)] <- range(q_grid[boot_res2[,1]<boot_res2[,3]])

TAB_a[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31),3] <- c(table_a[1:5,1], uptable_a[,1], table_a[6:11,1])
TAB_a[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32),3] <- c(table_a[1:5,2], uptable_a[,2], table_a[6:11,2])
TAB_a[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32),4] <- c(table_a[1:5,3], uptable_a[,3], table_a[6:11,3])
TAB_a[c(33,34),3] <- c(LT, CT)
TAB_a[35,3] <- gamma_hat
TAB_a[36,c(3,4)] <- range(q_grid[boot_res2[,1]<boot_res2[,3]])

CI<-range(q_grid[boot_res2[,1]<boot_res2[,3]])
ix1<-which(q_grid%in%CI)[1]
ix2<-which(q_grid%in%CI)[2]

pdf("grid_LEV.pdf", width=5, height=4.5)
par(mar = c(4,2,0.5,2))
plot(q_grid, boot_res2[,1], type="l", ylab="", cex.lab=1, xlab="Leverage")
lines(q_grid, boot_res2[,3], col="red", lty=2)
lines(x=c(CI[1],CI[1]), y=c(-10,boot_res2[ix1,3]), col="blue", lty=2)
lines(x=c(CI[2],CI[2]), y=c(-10,boot_res2[ix2,3]), col="blue", lty=2)
arrows(CI[1],0,CI[2],0, col="blue", length = 0.05)
arrows(CI[2],0,CI[1],0, col="blue", length = 0.05)
dev.off()

load("TQ_lin.Rdata")
load("TQ.Rdata")
cvs_boot_l <- apply((boot_res1),2,quantile,0.025)
cvs_boot_u <- apply((boot_res1),2,quantile,0.975)
table <- cbind(theta_hat, theta_hat-cvs_boot, theta_hat+cvs_boot)
table_a <- cbind(theta_hat, theta_hat-cvs_boot_u, theta_hat-cvs_boot_l)
upcvs_boot_l <- apply((boot_res1[,1:p1]+boot_res1[,(p1+2):p]),2,quantile,0.025)
upcvs_boot_u <- apply((boot_res1[,1:p1]+boot_res1[,(p1+2):p]),2,quantile,0.975)
alpha_up <- alpha_hat[1:p1]+alpha_hat[(p1+2):p]
uptable <- cbind(alpha_up, alpha_up-upcvs_boot, alpha_up+upcvs_boot)
uptable_a <- cbind(alpha_up, alpha_up-upcvs_boot_u, alpha_up-upcvs_boot_l)
TAB[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31),6] <- c(table[1:5,1], uptable[,1], table[6:11,1])
TAB[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32),6] <- c(table[1:5,2], uptable[,2], table[6:11,2])
TAB[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32),7] <- c(table[1:5,3], uptable[,3], table[6:11,3])
TAB[c(33,34),6] <- c(LT, CT)
TAB[35,6] <- gamma_hat
TAB[36,c(6,7)] <- range(q_grid[boot_res2[,1]<boot_res2[,3]])

TAB_a[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31),6] <- c(table_a[1:5,1], uptable_a[,1], table_a[6:11,1])
TAB_a[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32),6] <- c(table_a[1:5,2], uptable_a[,2], table_a[6:11,2])
TAB_a[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32),7] <- c(table_a[1:5,3], uptable_a[,3], table_a[6:11,3])
TAB_a[c(33,34),6] <- c(LT, CT)
TAB_a[35,6] <- gamma_hat
TAB_a[36,c(6,7)] <- range(q_grid[boot_res2[,1]<boot_res2[,3]])

CI<-range(q_grid[boot_res2[,1]<boot_res2[,3]])
ix1<-which(q_grid%in%CI)[1]
ix2<-which(q_grid%in%CI)[2]

pdf("grid_TQ2.pdf", width=5, height=3)
par(mar = c(4,2,0.5,2))
plot(q_grid, boot_res2[,1], type="l", ylab="", cex.lab=1, xlab="Tobin's Q")
lines(q_grid, boot_res2[,3], col="red", lty=2)
lines(x=c(CI[1],CI[1]), y=c(-10,boot_res2[ix1,3]), col="blue", lty=2)
lines(x=c(CI[2],CI[2]), y=c(-10,boot_res2[ix2,3]), col="blue", lty=2)
arrows(CI[1],0,CI[2],0, col="blue", length = 0.05)
arrows(CI[2],0,CI[1],0, col="blue", length = 0.05)
dev.off()

load("TQ_submin_lin.Rdata")
load("TQ_submin.Rdata")

cvs_boot_l <- apply((boot_res1),2,quantile,0.025)
cvs_boot_u <- apply((boot_res1),2,quantile,0.975)
table <- cbind(theta_hat, theta_hat-cvs_boot, theta_hat+cvs_boot)
table_a <- cbind(theta_hat, theta_hat-cvs_boot_u, theta_hat-cvs_boot_l)
upcvs_boot_l <- quantile((boot_res1[,p1]+boot_res1[,(p1+2):p]),0.025)
upcvs_boot_u <- quantile((boot_res1[,p1]+boot_res1[,(p1+2):p]),0.975)
upcvs_boot <- quantile(abs(boot_res1[,p1]+boot_res1[,(p1+2):p]),0.95)
alpha_up <- alpha_hat[p1]+alpha_hat[(p1+2):p]
uptable <- cbind(alpha_up, alpha_up-upcvs_boot, alpha_up+upcvs_boot)
uptable_a <- cbind(alpha_up, alpha_up-upcvs_boot_u, alpha_up-upcvs_boot_l)
TAB2[c(1,3,5,7,9,11,13,15),3] <- c(table[1:5,1], uptable[,1], table[6:7,1])
TAB2[c(2,4,6,8,10,12,14,16),3] <- c(table[1:5,2], uptable[,2], table[6:7,2])
TAB2[c(2,4,6,8,10,12,14,16),4] <- c(table[1:5,3], uptable[,3], table[6:7,3])
TAB2[c(17,18),3] <- c(LT, CT)
TAB2[19,3] <- gamma_hat
TAB2[20,c(3,4)] <- range(q_grid[boot_res2[,1]<boot_res2[,3]])

TAB2_a[c(1,3,5,7,9,11,13,15),3] <- c(table_a[1:5,1], uptable_a[,1], table_a[6:7,1])
TAB2_a[c(2,4,6,8,10,12,14,16),3] <- c(table_a[1:5,2], uptable_a[,2], table_a[6:7,2])
TAB2_a[c(2,4,6,8,10,12,14,16),4] <- c(table_a[1:5,3], uptable_a[,3], table_a[6:7,3])
TAB2_a[c(17,18),3] <- c(LT, CT)
TAB2_a[19,3] <- gamma_hat
TAB2_a[20,c(3,4)] <- range(q_grid[boot_res2[,1]<boot_res2[,3]])

CI<-range(q_grid[boot_res2[,1]<boot_res2[,3]])
ix1<-which(q_grid%in%CI)[1]
ix2<-which(q_grid%in%CI)[2]

pdf("grid_TQ1.pdf", width=5, height=3)
par(mar = c(4,2,0.5,2))
plot(q_grid, boot_res2[,1], type="l", ylab="", cex.lab=1, xlab="Tobin's Q")
lines(q_grid, boot_res2[,3], col="red", lty=2)
lines(x=c(CI[1],CI[1]), y=c(-10,boot_res2[ix1,3]), col="blue", lty=2)
lines(x=c(CI[2],CI[2]), y=c(-10,boot_res2[ix2,3]), col="blue", lty=2)
arrows(CI[1],0,CI[2],0, col="blue", length = 0.05)
arrows(CI[2],0,CI[1],0, col="blue", length = 0.05)
dev.off()



print(xtable(TAB_a, digits=3), type="latex", sanitize.text.function = function(x){x}, include.colnames=FALSE, include.rownames=FALSE)
print(xtable(TAB2_a, digits=3), type="latex", sanitize.text.function = function(x){x}, include.colnames=FALSE, include.rownames=FALSE)
