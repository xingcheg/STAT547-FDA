library(fdapace)
library(ggplot2)

############# functions ###############
## plot sparse (raw) functional data
plot.spar.fdata <- function(Lt, Ly){
  n <- length(Lt)
  x <- unlist(Lt)
  y <- unlist(Ly)
  count <- unlist(lapply(Lt, length))
  label <- as.factor(rep(1:n, times = count))
  DD <- data.frame(x = x, y = y, label = label)
  p <- ggplot(data = DD, aes(x = x, y = y, group = label)) + 
    geom_line(size = 0.5, alpha = 0.8) + 
    geom_point(colour = "blue", alpha = 0.25) + 
    theme_bw()
  print(p)
}

## evaluate 3 Methods for selecting K (for dense data).
den.model.eval.K <- function(fdata){
  input <- MakeFPCAInputs(tVec = fdata$pts, yVec = fdata$Y)
  fit1 <- FPCA(Ly = input$Ly, Lt = input$Lt,
               optns = list(methodSelectK = "FVE",
                            methodXi = "IN"))
  fit2 <- FPCA(Ly = input$Ly, Lt = input$Lt,
               optns = list(methodSelectK = "AIC",
                            methodXi = "IN"))
  fit3 <- FPCA(Ly = input$Ly, Lt = input$Lt,
               optns = list(methodSelectK = "BIC",
                            methodXi = "IN"))
  K1 <- fit1$selectK
  K2 <- fit2$selectK
  K3 <- fit3$selectK
  cat("FVE selection: K=", K1, "components; ",
      "AIC selection: K=", K2, "components; ",
      "BIC selection: K=", K3, "components.\n")
}

## evaluate 3 Methods for selecting K (for sparse data).
spar.model.eval.K <- function(fdata){
  fit1 <- FPCA(Ly = fdata$Ly, Lt = fdata$Lt,
               optns = list(methodSelectK = "FVE",
                            methodXi = "CE"))
  fit2 <- FPCA(Ly = fdata$Ly, Lt = fdata$Lt,
               optns = list(methodSelectK = "AIC",
                            methodXi = "CE"))
  fit3 <- FPCA(Ly = fdata$Ly, Lt = fdata$Lt,
               optns = list(methodSelectK = "BIC",
                            methodXi = "CE"))
  K1 <- fit1$selectK
  K2 <- fit2$selectK
  K3 <- fit3$selectK
  cat("FVE selection: K=", K1, "components; ",
      "AIC selection: K=", K2, "components; ",
      "BIC selection: K=", K3, "components.\n")
}

## evaluate 2 Methods for computing FPCs (for dense data)
den.model.eval.xi <- function(fdata){
  input <- MakeFPCAInputs(tVec = fdata$pts, yVec = fdata$Y)
  fit1 <- FPCA(Ly = input$Ly, Lt = input$Lt,
               optns = list(methodSelectK = "FVE",
                            methodXi = "IN"))
  fit2 <- FPCA(Ly = input$Ly, Lt = input$Lt,
               optns = list(methodSelectK = "FVE",
                            methodXi = "CE"))
  par(mfrow = c(2,3))
  plot(fit1$xiEst[,1], fit2$xiEst[,1], xlab = "Integration method",
       ylab = "BLUP", main = "1st FPC")
  plot(fit1$xiEst[,2], fit2$xiEst[,2], xlab = "Integration method",
       ylab = "BLUP", main = "2nd FPC")
  plot(fit1$xiEst[,3], fit2$xiEst[,3], xlab = "Integration method",
       ylab = "BLUP", main = "3rd FPC")
  matplot(fit1$workGrid, fit1$phi[,1:3], type = "l", xlab = "t", ylab = "First 3 Phi",
          main = "Integration method")
  matplot(fit2$workGrid, fit2$phi[,1:3], type = "l",  xlab = "t", ylab = "First 3 Phi",
          main = "BLUP")
  par(mfrow = c(1,1))
}

## evaluate 2 Methods for computing FPCs (for sparse data)
spar.model.eval.xi <- function(fdata){
  fit <- FPCA(Ly = fdata$Ly, Lt = fdata$Lt,
               optns = list(methodSelectK = "AIC",
                            methodXi = "CE"))
  
  par(mfrow = c(1,3))
  plot(fdata$xi[,1], fit$xiEst[,1], xlab = "TRUE",
       ylab = "BLUP", main = "1st FPC")
  plot(fdata$x[,2], fit$xiEst[,2], xlab = "TRUE",
       ylab = "BLUP", main = "2nd FPC")
  plot(fdata$x[,3], fit$xiEst[,3], xlab = "TRUE",
       ylab = "BLUP", main = "3rd FPC")
  par(mfrow = c(1,1))
}

########### check basis ##################
tt <- seq(0, 1, length.out = 101)
base3 <- CreateBasis(K = 3, pts = tt)
base20 <- CreateBasis(K = 20, pts = tt)
matplot(tt, base3, type = "l")
matplot(tt, base20, type = "l")

########### generate dense data ###############
## K = 3, sigma = 0.05
fdata_den_3_1 <- MakeGPFunctionalData(n = 100, M = 101, K = 3, sigma = 0.05,
                                      lambda = c(1, 0.8, 0.6))
## K = 3, sigma = 0.5
fdata_den_3_2 <- MakeGPFunctionalData(n = 100, M = 101, K = 3, sigma = 0.5,
                                      lambda = c(1, 0.8, 0.6))
## K = 20, sigma = 0.05
fdata_den_20_1 <- MakeGPFunctionalData(n = 100, M = 101, K = 20, sigma = 0.05,
                                      lambda = exp(-(1/2)*(0:19)) )
## K = 20, sigma = 0.5
fdata_den_20_2 <- MakeGPFunctionalData(n = 100, M = 101, K = 20, sigma = 0.5,
                                      lambda = exp(-(1/2)*(0:19)) )


matplot(tt, t(fdata_den_3_1$Y), type = "l")
matplot(tt, t(fdata_den_3_2$Y), type = "l")
matplot(tt, t(fdata_den_20_1$Y), type = "l")
matplot(tt, t(fdata_den_20_2$Y), type = "l")




########### generate sparse data ###############
## K = 3, sigma = 0.05
fdata_spar_3_1 <- MakeSparseGP(n = 100, sparsity = 6:12, K = 3,
                               lambda = c(1, 0.8, 0.6), sigma = 0.05)
## K = 3, sigma = 0.5
fdata_spar_3_2 <- MakeSparseGP(n = 100, sparsity = 6:12, K = 3,
                               lambda = c(1, 0.8, 0.6), sigma = 0.5)
## K = 20, sigma = 0.05
fdata_spar_20_1 <- MakeSparseGP(n = 100, sparsity = 6:12, K = 20,
                               lambda =  exp(-(1/2)*(0:19)), sigma = 0.05)
## K = 20, sigma = 0.5
fdata_spar_20_2 <- MakeSparseGP(n = 100, sparsity = 6:12, K = 20,
                               lambda =  exp(-(1/2)*(0:19)), sigma = 0.5)



plot.spar.fdata(fdata_spar_3_1$Lt, fdata_spar_3_1$Ly)
plot.spar.fdata(fdata_spar_3_2$Lt, fdata_spar_3_2$Ly)

plot.spar.fdata(fdata_spar_20_1$Lt, fdata_spar_20_1$Ly)
plot.spar.fdata(fdata_spar_20_2$Lt, fdata_spar_20_2$Ly)



################ FPCA for dense functional data ##################
input_den_3_1 <- MakeFPCAInputs(tVec = tt, yVec = fdata_den_3_1$Y)
fit_den_3_1 <- FPCA(Ly = input_den_3_1$Ly, Lt = input_den_3_1$Lt,
                    optns = list(methodSelectK = 'FVE',
                                 methodXi = "IN"))



###### evaluate K selection methods ###########

## dense, K = 3, sigma = 0.05
den.model.eval.K(fdata_den_3_1)
## dense, K = 3, sigma = 0.5
den.model.eval.K(fdata_den_3_2)
## dense, K = 20, sigma = 0.05
den.model.eval.K(fdata_den_20_1)
## dense, K = 20, sigma = 0.5
den.model.eval.K(fdata_den_20_2)

## sparse, K = 3, sigma = 0.05
spar.model.eval.K(fdata_spar_3_1)
## sparse, K = 3, sigma = 0.5
spar.model.eval.K(fdata_spar_3_2)
## sparse, K = 20, sigma = 0.05
spar.model.eval.K(fdata_spar_20_1)
## sparse, K = 20, sigma = 0.5
spar.model.eval.K(fdata_spar_20_2)





###### integration method vs BLUP ###########

## dense, K = 3, sigma = 0.05
den.model.eval.xi(fdata_den_3_1)
## dense, K = 3, sigma = 0.5
den.model.eval.xi(fdata_den_3_2)
## dense, K = 20, sigma = 0.05
den.model.eval.xi(fdata_den_20_1)
## dense, K = 20, sigma = 0.5
den.model.eval.xi(fdata_den_20_2)

## sparse, K = 3, sigma = 0.05
spar.model.eval.xi(fdata_spar_3_1)
## sparse, K = 3, sigma = 0.5
spar.model.eval.xi(fdata_spar_3_2)
## sparse, K = 20, sigma = 0.05
spar.model.eval.xi(fdata_spar_20_1)
## sparse, K = 20, sigma = 0.5
spar.model.eval.xi(fdata_spar_20_2)



