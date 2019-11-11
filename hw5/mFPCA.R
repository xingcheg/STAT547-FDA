# Note: These functions are used to estimate bivariate FPCA (sparse).
# The methodology applied here are based on paper: Chiou, Chen and Yang (2014). 
# I also import some functions in R package: fdapace, locfit.
# The parameters (bindwidths, K) can only be manually adjusted. 

library(fdapace)
library(locfit) 
library(dplyr)

## create mFPCA object that can be directed used in mFPCA function.
Make_mFPCA_Inputs <- function(idVec, tVec, yMat){
  label <- unique(idVec)
  n_sub <- length(label)
  count <- rep(0, n_sub)
  for (i in 1:n_sub){
    count[i] <- sum(idVec == label[i])
  }
  weight <- rep(1/count, times = count)
  out <- data.frame(id = idVec, t = tVec, weight, yMat)
  return(out)
}


## estimate mean funcitons (D: output of Make_mFPCA_Inputs).
Mean_Curve_Est <- function(D, deg = 1, h_Mu = 0.1, ngrid = 100){
  N <- nrow(D)
  p <- ncol(D) - 3
  t <- D$t
  Y <- D[,-(1:3)]
  muhat_grid <- matrix(0, ngrid, p)
  muhat <- matrix(0, N, p)
  for (i in 1:p){
    fit <- locfit(Y[,i] ~ lp(t, deg = deg, h = h_Mu), weights = D$weight,
                  ev = lfgrid(mg = ngrid) )
    muhat_grid[,i] <- predict(fit)
    muhat[,i] <- predict(fit, t)
  }
  return(list(muhat = muhat, muhat_grid = muhat_grid))
}


## create raw covariance data (D: output of Make_mFPCA_Inputs).
Raw_Cov_data <- function(D, muhat){
  D[,-(1:3)] <-  D[,-(1:3)] - muhat
  p <- ncol(muhat)
  label <- unique(D$id)
  n_sub <- length(label)
  rcov <- NULL
  for (i in 1:n_sub){
    d <- D[ D$id==label[i], ]
    TT <- expand.grid(t1=d$t, t2=d$t) 
    w <- d$weight[1]
    if ( w ==1 ){
      weight <- 1/2
    } else {
      weight <- rep( (w^2)/(1-w), nrow(TT) )
    }
    id <- rep(label[i], nrow(TT))
    Raw <- apply(d[,-(1:3)], 2, function(x){
      c(tcrossprod(x))
    })
    if (nrow(TT) == 1){
      dd <- c(id, as.numeric(TT), weight, as.numeric(Raw))
    } else {
      dd <- cbind(id, TT, weight, Raw)
    }
    rcov <- rbind(rcov, dd)
  }
  rcov <- as.data.frame(rcov)
  return(rcov)
}


## estimate covariance surface (rcov: output of Raw_Cov_data).
Cov_Surf_Est <- function(rcov, t, h_Cov = 0.2, ngrid = 100){
  p <- ncol(rcov) - 4
  rcov <- rcov %>% filter(t1 != t2)
  t1 <- rcov$t1
  t2 <- rcov$t2
  weight <- rcov$weight
  covArray <- array(0, c(ngrid, ngrid, p))
  covArray_raw <- array(0, c(ngrid, ngrid, p))
  covDiag <- matrix(0, length(t), p)
  covDiag_grid <- matrix(0, ngrid, p)
  for (i in 1:p){
    yi <- rcov[,4+i]
    fit <- locfit(yi ~ lp(t1, t2, deg=1, h=h_Cov), 
                  ev=lfgrid(mg=ngrid), weights = weight)
    covMat <- matrix(predict(fit), ngrid, ngrid) 
    eig_cov <- eigen(covMat)
    eig_value <- eig_cov$values
    eig_vector <- eig_cov$vectors
    eff_ind <- which(eig_value > 0)
    covMat <- eig_vector[,eff_ind] %*%  (t(eig_vector[,eff_ind]) *  eig_value[eff_ind])
    covDiag[,i] <- pmax( predict(fit, data.frame(t1=t, t2=t)), 1e-4)
    covDiag_grid[,i] <- diag(covMat)
    V <- diag(diag(covMat)^(-1/2))
    covArray[,,i] <- V %*% covMat %*% V
    covArray_raw[,,i] <- covMat
  }
  return(list(covArray = covArray, covDiag = covDiag, covDiag_grid = covDiag_grid,
              covArray_raw = covArray_raw))
}


## data transformation (D: output of Make_mFPCA_Inputs).
Trans_data <- function(D, muhat, covDiag){
  D[,-(1:3)] <-  (D[,-(1:3)] - muhat) / ((covDiag)^(1/2))
  return(D)
}


## create raw cross-covariance data (D: output of Trans_data).
Raw_Cross_Cov_data <- function(D){
  p <- ncol(D) - 3
  label <- unique(D$id)
  n_sub <- length(label)
  rcrosscov <- NULL
  for (i in 1:n_sub){
    d <- D[ D$id==label[i], ]
    TT <- expand.grid(t1=d$t, t2=d$t) 
    w <- d$weight[1]
    if ( w ==1 ){
      weight <- 1/2
    } else {
      weight <- rep( (w^2)/(1-w), nrow(TT) )
    }
    id <- rep(label[i], nrow(TT))
    Raw <- matrix(0, nrow(TT), p*(p-1)/2 )
    flag <- 0
    for (k in 2:p){
      for (l in 1:(k-1)){
        flag <- flag + 1
        Raw[,flag] <- c( d[,3+k] %*% t(d[,3+l]) )
      }
    }
    if (nrow(TT) == 1){
      dd <- c(id, as.numeric(TT), weight, as.numeric(Raw))
    } else {
      dd <- cbind(id, TT, weight, Raw)
    }
    rcrosscov <- rbind(rcrosscov, dd)
  }
  rcrosscov <- as.data.frame(rcrosscov)
  return(rcrosscov)
}


## estimate cross-covariance surface (rcrosscov: output of Raw_Cross_Cov_data).
Cross_Cov_Surf_Est <- function(rcrosscov, h_Cov = 0.2, ngrid = 100){
  q <- ncol(rcrosscov) - 4 
  rcrosscov <- rcrosscov %>% filter(t1 != t2)
  t1 <- rcrosscov$t1
  t2 <- rcrosscov$t2
  weight <- rcrosscov$weight
  crosscovArray <- array( 0, c(ngrid, ngrid, q) )
  flag <- 0
  for (i in 1:q){
      flag <- flag + 1
      y_kl <- rcrosscov[,4+i]
      fit <- locfit(y_kl ~ lp(t1, t2, deg=1, h=h_Cov), 
                    ev=lfgrid(mg=ngrid), weights = weight)
      crosscovMat <- matrix(predict(fit), ngrid, ngrid) 
      crosscovArray[,,i] <- crosscovMat
    }

  return(crosscovArray)
}


## estimate sigma_k^2 (D: output of Make_mFPCA_Inputs).
Sigma_Est <- function(D, muhat, covArray_raw, h_W = 0.2){
  D[,-(1:3)] <-  D[,-(1:3)] - muhat
  p <- ncol(D) - 3
  tt <- D$t
  weight <- D$weight
  ngrid <- dim(covArray_raw)[1]
  eff_index <- seq( round(ngrid/4), round(3*ngrid/4), 1 )
  sigma2 <- rep(0, p) 
  for (i in 1:p){
    yi <- (D[,3+i])^2
    fit <- locfit(yi ~ lp(tt, deg=1, h=h_W), 
                  ev=lfgrid(mg=ngrid), weights = weight)
    errVec <- predict(fit)
    covVec <- diag(covArray_raw[,,i])
    sigmaVec <- errVec - covVec
    sigma2[i] <- mean( sigmaVec[eff_index] )
  }
  return(sigma2)
}


## compute first L eigenvalues and eigenfunctions by discretizing the
## covariance and cross-covariance functions
Get_Eigen <- function(L = 3, covArray, crosscovArray){
  p <- dim(covArray)[3]
  ngrid <- dim(covArray)[1]
  Cov_joint <- matrix(0, p*ngrid, p*ngrid)
  Cov_joint[1:ngrid, 1:ngrid] <- covArray[,,1]
  flag <- 0
  for (k in 2:p){
    Cov_joint[ngrid*(k-1)+(1:ngrid), ngrid*(k-1)+(1:ngrid)] <- covArray[,,k]
    for (l in 1:(k-1)){
      flag <- flag + 1
      Cov_joint[ngrid*(k-1)+(1:ngrid), ngrid*(l-1)+(1:ngrid)] <- crosscovArray[,,flag]
      Cov_joint[ngrid*(l-1)+(1:ngrid), ngrid*(k-1)+(1:ngrid)] <- t(crosscovArray[,,flag])
    }
  }
  eigenCov <- eigen(Cov_joint)
  L1 <- which(eigenCov$values<0)[1]
  L <- min(L, L1-1)
  eigenval <- eigenCov$values[1:L]
  eigenfunc <- eigenCov$vectors[,1:L, drop = FALSE]
  Cov_joint_pd <- eigenfunc %*% diag(eigenval, ncol = length(eigenval)) %*% t(eigenfunc)
  return(list(eigenval = eigenval, eigenfunc = eigenfunc,
              Cov_joint = Cov_joint,  Cov_joint_pd = Cov_joint_pd))
}


## BLUP (D: output of Trans_data, eigen_out: output of Get_Eigen).
mBLUP <- function(D, eigen_out, sigma2, covDiag, tgrid){
  p <- ncol(D) - 3
  L <- length(eigen_out$eigenval)
  label <- unique(D$id)
  n_sub <- length(label)
  Xi <- matrix(0, n_sub, L)
  rownames(Xi) <- label
  ngrid <- length(tgrid)
  for (i in 1:n_sub){
    index_i <- which(D$id==label[i])
    d <- D[ index_i, ]
    covDiag_i <- covDiag[ index_i, , drop = FALSE]
    Ui <- c( as.matrix(d[,-(1:3)]) )
    ti <- d$t
    mi <- length(ti)
    Hi_trans <- matrix(0, p*mi, L)
    Sig_Ui <- matrix(0, p*mi, p*mi)
    
    Hi_trans[1:mi, ] <- matrix(ConvertSupport(tgrid, ti, 
      phi = eigen_out$eigenfunc[1:ngrid, ]), ncol = L)
    Sig_Ui[1:mi, 1:mi] <- ConvertSupport(tgrid, ti,
      Cov = eigen_out$Cov_joint_pd[1:ngrid, 1:ngrid]) + 
      diag(sigma2[1]/covDiag_i[,1], nrow = mi)
    
    for (k in 2:p){
      Hi_trans[(k-1)*mi+(1:mi),] <- matrix(ConvertSupport(tgrid, ti, 
        phi = eigen_out$eigenfunc[(k-1)*ngrid+(1:ngrid),]), ncol = L)
      Sig_Ui[(k-1)*mi+(1:mi), (k-1)*mi+(1:mi)] <- ConvertSupport(tgrid, ti,
        Cov = eigen_out$Cov_joint_pd[(k-1)*ngrid+(1:ngrid), (k-1)*ngrid+(1:ngrid)]) + 
        diag( sigma2[k] / covDiag_i[, k], nrow = mi)
      for (l in 1:(k-1)){
        Sig_Ui[(k-1)*mi+(1:mi), (l-1)*mi+(1:mi)] <- ConvertSupport(tgrid, ti,
          Cov = eigen_out$Cov_joint_pd[(k-1)*ngrid+(1:ngrid), (l-1)*ngrid+(1:ngrid)], isCrossCov = TRUE)
        Sig_Ui[(l-1)*mi+(1:mi), (k-1)*mi+(1:mi)] <- t(Sig_Ui[(k-1)*mi+(1:mi), (l-1)*mi+(1:mi)])
      }
    }
    Hi <- t(Hi_trans) * eigen_out$eigenval
    Xi[i,] <- as.numeric( Hi %*% solve(Sig_Ui) %*% Ui )
  }
  return(Xi)
}


## estimate X_L
X_L_est <- function(muhat_grid, eigen_out, XiEst, covDiag_grid){
  ngrid <- nrow(muhat_grid)
  p <- ncol(muhat_grid)
  n_sub <- nrow(XiEst)
  Xhat <- array(0, c(n_sub, ngrid, p))
  for (i in 1:n_sub){
    Xhat[i,,] <- muhat_grid
  }
  Phi <- eigen_out$eigenfunc
  for (k in 1:p){
    Phi_k <- Phi[(k-1)*ngrid+(1:ngrid),]
    D_Phi_k <- Phi_k * sqrt(covDiag_grid[,k])
    Xhat[,,k] <- Xhat[,,k] + matrix( XiEst %*% t(D_Phi_k), ncol = ngrid)
  }
  return(Xhat)
}


## mFPCA
 mFPCA <- function(idVec, tVec, yMat, ngrid = 100, h_Mu = 0.1, h_Cov = 0.2,
                   h_W = 0.2, L = 3){
   DD1 <- Make_mFPCA_Inputs(idVec, tVec, yMat)
   tgrid <- seq(min(DD1$t), max(DD1$t), length.out = ngrid)
   
   ## estimate mean curves
   muEst <- Mean_Curve_Est(DD1, h_Mu = h_Mu)
   muhat <- muEst$muhat
   muhat_grid <- muEst$muhat_grid
   
   ## estimate covariance surfaces
   rcov <- Raw_Cov_data(DD1, muhat)
   covEst <- Cov_Surf_Est(rcov, DD1$t, h_Cov = h_Cov)
   covArray <- covEst$covArray
   covDiag <- covEst$covDiag
   covDiag_grid <- covEst$covDiag_grid
   covArray_raw <- covEst$covArray_raw
   
   ## transform data & estimate cross-covariance surfaces 
   DD2 <- Trans_data(DD1, muhat, covDiag)
   rcrosscov <- Raw_Cross_Cov_data(DD2)
   crosscovArray <- Cross_Cov_Surf_Est(rcrosscov, h_Cov = h_Cov)
   sigma2 <- Sigma_Est(DD1, muhat, covArray_raw, h_W = h_W)
   
   ## eigenvalue & eigenfunction 
   eigen_out <- Get_Eigen(L = L, covArray, crosscovArray)
   Cov_joint <- eigen_out$Cov_joint
   eigenfunc <- eigen_out$eigenfunc
   
   ## BLUP for estimating xi
   XiEst <- mBLUP(DD2, eigen_out, sigma2, covDiag, tgrid)
   
   ## Estimate X_L 
   Xhat <- X_L_est(muhat_grid, eigen_out, XiEst, covDiag_grid)
   
   output <- list(Mu_hat = muhat_grid, D_hat = covDiag_grid, C_hat = eigen_out$Cov_joint_pd,
                  Phi_hat = eigenfunc, Xi_hat = XiEst, sigma2_hat = sigma2, X_hat = Xhat, tgrid = tgrid)
   
   return(output)
 }



