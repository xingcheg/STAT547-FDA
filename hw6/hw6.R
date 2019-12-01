setwd("/Users/apple/Desktop/ISU 2019 fall/STAT547/hw/hw6")
library(fdapace)
library(ggplot2)
library(tidyverse)
library(fda)

######################### Problem 1 #######################
## scalar-response functional linear regression model
###########################################################
load("medfly25.RData")
source("sca_func_reg.R")
dat <- medfly25 

## visualization
ggplot(data = dat) + 
  geom_line(aes(x = Days, y = nEggs, group = ID), size = 0.4, alpha = 0.6) + 
  theme_bw()

ggplot(data = dat[dat$ID %in% 1:5,]) + 
  geom_line(aes(x = Days, y = nEggs, colour = as.factor(ID))) + 
  theme_bw()


## X, Y data
tGrid <- 1:25
denseSamp <- MakeFPCAInputs(IDs=dat$ID, tVec=dat$Days, yVec=dat$nEggs) 
X <- matrix( as.numeric(unlist(denseSamp$Ly)), ncol = length(denseSamp$Ly) )
Y <- dat %>%  
  group_by(ID) %>%  
  summarize(y=nEggsRemain[1]) %>%  
  .$y 


## regularization method 
res_out <- sca_func_reg_pen(X, Y, 1:25, nbasis = 7, norder = 4, 
                            dropind = 1, pen_order = 2)

## FPC regression 
res_out1 <- FPC_reg(X = denseSamp, Y, 1:25, optnsX=list(), K = 10)


plot(res_out$Yhat, res_out1$yhat, xlab = "Regularization Method", ylab = "FPC Method")


rm(list=ls())





######################### Problem 2 #######################
## functional concurrent regression model
###########################################################

dat <- read.csv("USGDP.csv")
totalPopulation <- as.double( dat$totalPopulation )
totalLaborForce <- as.double( dat$totalLaborForce )
perCapitaGDP <- as.double( dat$perCapitaGDP )



## Z1, Z2, Y data
tGrid <- 1997:2015
n <- 49
Z1 <-  matrix(totalPopulation, ncol = n, byrow = TRUE) 
Z2 <-  matrix(totalLaborForce, ncol = n, byrow = TRUE) 
Y <-  matrix(perCapitaGDP, ncol = n, byrow = TRUE) 


## X1 <- log(Z1); X2 <- Z2 / Z1; Y <- log(Y)
X1 <- log(Z1)
X2 <- Z2 / Z1
Y <- log(Y)

par(mfrow = c(2,2))
matplot(tGrid, X1, type = "l")
matplot(tGrid, X2, type = "l")
matplot(tGrid, Y, type = "l")
par(mfrow = c(1,1))

## FCReg
X1_list <- MakeFPCAInputs(tVec = tGrid, yVec = t(X1)) 
X2_list <- MakeFPCAInputs(tVec = tGrid, yVec = t(X2)) 
Y_list <- MakeFPCAInputs(tVec = tGrid, yVec = t(Y)) 

FCR_res <- FCReg(vars = list(X1 = X1_list, X2 = X2_list, Y = Y_list), userBwMu = 3, 
                 userBwCov = 6, outGrid = tGrid)


par(mfrow = c(2,2))
plot(tGrid, FCR_res$beta0, type = "l", ylab = "beta0", main = "beta0")
plot(tGrid, FCR_res$beta[1,], type = "l", ylab = "beta1", main = "beta1")
plot(tGrid, FCR_res$beta[2,], type = "l", ylab = "beta2", main = "beta2")

matplot(tGrid, FCR_res$beta0 %*% t(rep(1, n)) + 
          FCR_res$beta[1,] * X1 + 
          FCR_res$beta[2,] * X2, type = "l", ylab = "yhat", main = "yhat")

par(mfrow = c(1,1))





######################### Problem 3 #######################
## function-on-function regression model (FPC method)
###########################################################

X1_list1 <- list(Lt = X1_list$Lt, Ly = X1_list$Ly)
X2_list1 <- list(Lt = X2_list$Lt, Ly = X2_list$Ly)
Y_list1 <- list(Lt = Y_list$Lt, Ly = Y_list$Ly)

FPC_res <- FPCReg(vars = list(X1 = X1_list1, X2 = X2_list1, Y = Y_list1))
tGrid_fine <- seq(1997, 2015, length.out = 51)

par(mfrow = c(2,2))
beta <- FPC_res$estiBeta
image(tGrid_fine, tGrid_fine, beta$betaX1Y, main = "beta1")
image(tGrid_fine, tGrid_fine, beta$betaX2Y, main = "beta2")

yhat <- matrix(unlist(FPC_res$predictY), 51)
matplot(tGrid_fine, yhat, type = "l", main = "yhat")
par(mfrow = c(1,1))





