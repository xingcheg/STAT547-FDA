library(fda)
library(fdapace)

############### Problem 1 #################
data(gait)
Hip <- gait[, , 1]
Knee <- gait[, , 2]
time <- as.numeric( rownames(Hip) )
n <- ncol(Hip)
m <- nrow(Hip)

XX <- t(rbind(Hip, Knee))
mu <- colMeans(XX) 
GG <- cov(XX)

## plot the mean function
par(mfrow = c(1,2))
matplot(time, Hip, type='l', col='black', main = "Hip Angle", xlab = "time") 
lines(time, mu[1:m], col = "red", lwd = 2.5)

matplot(time, Knee, type='l', col='black', main = "Knee Angle", xlab = "time") 
lines(time, mu[(m+1):(2*m)], col = "red", lwd = 2.5)
par(mfrow = c(1,1))

## check the total variation
eig <- eigen(GG)
lam <- eig$values / m 
FVEeach <- lam / sum(lam) 
FVE <- cumsum(FVEeach) 

plot(FVE, pch = 16)
abline(h = 0.95, col = "red")

## plot the eigen-functions
phi <- eig$vectors * sqrt(m)
phi1 <- phi[1:m,]
phi2 <- phi[(m+1):(2*m),]

par(mfrow = c(1,2))
matplot(time, phi1[,1:2], type='l', col=c('red', 'blue'), xlab = "time", ylab = "hip") 
matplot(time, phi2[,1:2], type='l', col=c('red', 'blue'), xlab = "time", ylab = "knee") 
par(mfrow = c(1,1))

## Mode of variation plot 
par(mfcol = c(2,2))
lamphi11 <- sqrt(lam[1]) * phi1[, 1] 
plot(time, mu[1:m], type='l', lwd=2, ylim = c(-10, 60), 
     main = "FPC1", xlab = "time", ylab = "hip") 
lines(time, mu[1:m] + lamphi11, type='l', lty=2, col='red') 
lines(time, mu[1:m] - lamphi11, type='l', lty=2, col='blue') 

lamphi12 <- sqrt(lam[2]) * phi1[, 2] 
plot(time, mu[1:m], type='l', lwd=2, ylim = c(-10, 60), 
     main = "FPC2", xlab = "time", ylab = "hip") 
lines(time, mu[1:m] + lamphi12, type='l', lty=2, col='red') 
lines(time, mu[1:m] - lamphi12, type='l', lty=2, col='blue') 

lamphi21 <- sqrt(lam[1]) * phi2[, 1] 
plot(time, mu[(m+1):(2*m)], type='l', lwd=2, ylim = c(-0, 80), 
     main = "FPC1", xlab = "time", ylab = "knee") 
lines(time, mu[(m+1):(2*m)] + lamphi21, type='l', lty=2, col='red') 
lines(time, mu[(m+1):(2*m)] - lamphi21, type='l', lty=2, col='blue') 

lamphi22 <- sqrt(lam[2]) * phi2[, 2] 
plot(time, mu[(m+1):(2*m)], type='l', lwd=2, ylim = c(-0, 80), 
     main = "FPC2", xlab = "time", ylab = "knee") 
lines(time, mu[(m+1):(2*m)] + lamphi22, type='l', lty=2, col='red') 
lines(time, mu[(m+1):(2*m)] - lamphi22, type='l', lty=2, col='blue') 
par(mfrow = c(1,1))


## FPCs
XXcenter <- XX - matrix( rep(mu, n), nrow=n, ncol=2*m, byrow=TRUE) 
xi <- XXcenter %*% phi / m 

plot(xi[, 1], xi[, 2], xlab = "FPC1", ylab = "FPC2", 
     col = "red", cex = 0.001) 
text(xi[, 1], xi[, 2], label= sub("boy", "", rownames(XX))) 






############### Problem 2 #################
K <- 20 
n <- 50 
m <- 101 
t <- seq(0, 1, length.out = m) 
mu <- -(t - 0.5)^2 + 1 
plot(t, mu, type = "l") 



## Create Basis 
phiSim <- CreateBasis(K=K, t, type='fourier') 
matplot(t, phiSim[, 1:5], type='l') 
lamSim <- exp(-(1:K)) 

## Gaussian Process
# rows are for individuals, col are for different FPC 
xi <- matrix(rnorm(n * K), nrow = n, ncol = K) %*% diag(sqrt(lamSim)) 
X1 <- matrix(mu, nrow=n, ncol=m, byrow=TRUE) + xi %*% t(phiSim) 
matplot(t, t(X1), type='l') 


## non-Gaussian Process
xi <- matrix(rt(n * K, df = 2), nrow = n, ncol = K) %*% diag(sqrt(lamSim)) 
X2 <- matrix(mu, nrow=n, ncol=m, byrow=TRUE) + xi %*% t(phiSim) 
matplot(t, t(X2), type='l') 





############### Problem 3 #################
yeast <- read.table(file = "/Users/apple/Desktop/ISU 2019 fall/STAT547/data/yeast.txt")
time <- as.numeric( sub("alpha", "", names(yeast)) )

## original data
X <- as.matrix(yeast)
matplot(time, t(X), type='l', col=c(rep(2,44), rep(4,45)), 
        xlab = "time", ylab = "yeast") 

n <- nrow(X)
m <- ncol(X)

## smoothed data
X_smooth <- matrix(0, ncol = m, nrow = n)

for (i in 1:n){
  ind <- which(is.na(X[i,]))
  if (length(ind) > 0){
    fit <- smooth.spline(time[-ind], X[i,][-ind])
  }else{
    fit <- smooth.spline(time, X[i,], all.knots = TRUE, spar = 0.1)
  }
    xhat <- predict(fit, as.data.frame(time))$y$time
    X_smooth[i,] <- xhat
}

matplot(time, t(X_smooth), type='l', col=c(rep(2,44), rep(4,45)), 
        xlab = "time", ylab = "yeast") 


## check the total variation
mu <- colMeans(X_smooth) 
G <- cov(X_smooth) 
eig <- eigen(G) 
lam <- eig$values * diff(range(time)) / m 
FVEeach <- lam / sum(lam) 
FVE <- cumsum(FVEeach) 

plot(FVE, pch = 16)
abline(h = 0.95, col = "red")


## Mode of variation plot 
phi <- eig$vectors / sqrt(diff(range(time) / m)) 

par(mfrow = c(2,2))
lamphi1 <- sqrt(lam[1]) * phi[, 1] 
plot(time, mu[1:m], type='l', lwd=2, ylim = c(-1,1),
     main = "FPC1", xlab = "time", ylab = "PC1") 
lines(time, mu + lamphi1, type='l', lty=2, col='red') 
lines(time, mu - lamphi1, type='l', lty=2, col='blue') 

lamphi2 <- sqrt(lam[2]) * phi[, 2] 
plot(time, mu, type='l', lwd=2, ylim = c(-1, 1), 
     main = "FPC2", xlab = "time", ylab = "PC2") 
lines(time, mu + lamphi2, type='l', lty=2, col='red') 
lines(time, mu - lamphi2, type='l', lty=2, col='blue') 

lamphi3 <- sqrt(lam[3]) * phi[, 3] 
plot(time, mu, type='l', lwd=2, ylim = c(-1, 1), 
     main = "FPC3", xlab = "time", ylab = "PC3") 
lines(time, mu + lamphi3, type='l', lty=2, col='red') 
lines(time, mu - lamphi3, type='l', lty=2, col='blue') 

lamphi4 <- sqrt(lam[4]) * phi[, 4] 
plot(time, mu, type='l', lwd=2, ylim = c(-1, 1), 
     main = "FPC4", xlab = "time", ylab = "PC4") 
lines(time, mu + lamphi4, type='l', lty=2, col='red') 
lines(time, mu - lamphi4, type='l', lty=2, col='blue') 
par(mfrow = c(1,1))


## FPC
Xcenter <- X - matrix(mu, nrow=n, ncol=m, byrow=TRUE) 
xi <- Xcenter %*% phi * diff(range(time)) / m 
plot(xi[, 1], xi[, 2],  col=c(rep(2,44), rep(4,45)), pch = 16,
     xlab = "FPC1", ylab = "FPC2") 


