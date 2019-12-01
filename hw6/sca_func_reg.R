library(fda)


###########################################################
## scalar on functional regression (regularization)
###########################################################

sca_func_reg_pen <- function(X, Y, tGrid, nbasis, norder, dropind = NULL, pen_order = 2){
  par(mfrow = c(2,2))
  n <- length(Y)
  if (is.null(dropind)){
    bbasis <- create.bspline.basis(rangeval = range(tGrid), nbasis = nbasis, 
                                   norder = norder)
  } else {
    bbasis <- create.bspline.basis(rangeval = range(tGrid), nbasis = nbasis, 
                                   norder = norder, dropind = dropind)
  }
  bbasis_val <- eval.basis(tGrid, bbasis)
  X_fd <- smooth.basis(tGrid, X, bbasis)$fd
  plot(X_fd, main = "Smoothed X")
  
  C <- t(X_fd$coefs)
  if (is.null(dropind)){
    Omega <- bsplinepen(bbasis, Lfdobj = pen_order)
    J <- bsplinepen(bbasis, Lfdobj = 0)
  } else{
    Omega <- bsplinepen(bbasis, Lfdobj = pen_order)[-dropind, -dropind]
    J <- bsplinepen(bbasis, Lfdobj = 0)[-dropind, -dropind]
  }
  
  S <- C %*% J
  lambda <- 10^seq(-6, 6)
  n_lam <- length(lambda)
  gcv <- rep(0, n_lam)
  
  A12 <- apply( S, 2, sum )
  B <- rbind(rep(1, n), t(S))
  
  for (i in 1:n_lam){
    lam <- lambda[i]
    A22 <- crossprod( S ) + lam * Omega
    A <- rbind( c(n, A12), cbind(A12, A22) )
    H <- t(B) %*% solve(A) %*% B
    Yhat <- H %*% Y
    gcv[i] <- mean( (Y - Yhat)^2 ) / ( 1 - mean( diag(H) ) )^2
  }
  
  plot(log10(lambda), gcv, type = "b", main = "GCV", ylab = "gcv", xlab = "log10_lambda")
  lam1 <- lambda[which.min(gcv)]
  A22 <- crossprod( S ) + lam1 * Omega
  A <- rbind( c(n, A12), cbind(A12, A22) )
  est <- solve(A) %*% B %*% Y
  alpha <- est[1]
  b <- est[-1]
  beta <- bbasis_val %*% b
  Yhat <- t(B) %*% est
  
  plot(tGrid, beta, main = "Estimate of beta", type = "l")
  abline(h = 0, col = "red")
  plot(Y, Yhat, main = "Y vs Yhat")
  par(mfrow = c(1,1))
  
  return(list(Yhat = Yhat, alpha = alpha, beta = beta, lambda = lam1))
  
}





###########################################################
## scalar on functional regression (FPC regression)
###########################################################

FPC_reg <- function(X, Y, tGrid, optnsX=list(), K = NULL) { 
  
  par(mfrow = c(2,2))
  fpcaRes <- FPCA(X$Ly, X$Lt, optnsX) 
  if (is.null(K)){
    K <- ncol(fpcaRes$xiEst)
  } else {
    K <- min(K, ncol(fpcaRes$xiEst))
  }
  xiEst <- fpcaRes$xiEst[,1:K]
  dat <- data.frame(y=Y, x=xiEst) 
  model <- lm(y ~ ., dat) 
  b <- model$coefficients[-1] 
  a <- model$coefficients[1] 
  betaFun <- fpcaRes$phi[,1:K] %*% matrix(b) 
  alpha <- a 
  Yhat <- model$fitted.values
  Xhat <- fpcaRes$mu + fpcaRes$phi[,1:K] %*% t(fpcaRes$xiEst[,1:K])
  
  matplot(tGrid, Xhat, main = "Smoothed X", type = "l")
  plot(fpcaRes$cumFVE, main = "FVE", xlab = "K", ylab = "cumFVE")
  plot(tGrid, betaFun,  main = "Estimate of beta", type = "l")
  abline(h = 0, col = "red")
  plot(Y, Yhat, main = "Y vs Yhat")
  par(mfrow = c(1,1))
  
  res <- list( 
    alpha = alpha, 
    beta = betaFun,  
    xGrid = fpcaRes$workGrid, 
    yhat = Yhat
  ) 
  
  return(res)
} 
  
  
  
  
  
  
  
  
  
  









  
  
  
  
  


