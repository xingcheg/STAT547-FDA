library(dplyr) 
library(ggplot2) 
library(locfit) 
library(plotly)
library(fda)
library(reshape2)
source("/Users/apple/Desktop/ISU 2019 fall/STAT547/hw/hw5/mFPCA.R")

#### read data ####
DD0 <- read.delim('http://www.statsci.org/data/oz/wallaby.txt') %>%  
  select(Anim, Leng, Weight, Age) %>%
  na.omit %>% 
  # filter(Age >= 0.6 * 365.24 & Age <= 2.4 * 365.24 ) %>% 
  mutate(Leng = log(Leng), Weight = log(Weight)) %>%
  mutate(Leng = (Leng - mean(Leng))/sd(Leng), 
         Weight = (Weight - mean(Weight))/sd(Weight),
         Age = (Age - min(Age))/diff(range(Age)))

names(DD0) <- c("label", "y1", "y2", "t")


#### run mFPCA ####
idVec <- DD0$label
tVec <- DD0$t
yMat <- DD0[,c(2,3)]

mFPCA_out <- mFPCA(idVec, tVec, yMat, ngrid = 100, 
                   h_Mu = 0.2, h_Cov = 0.5, h_W = 0.5, L = 3)


#### plot mean ####
tgrid <- mFPCA_out$tgrid
muhat_grid <- mFPCA_out$Mu_hat
DD_mu <- data.frame(tgrid, muhat_grid)
names(DD_mu) <- c("tgrid", "mu1", "mu2")

# y1
ggplot(data = DD0) + 
  geom_line( aes(x = t, y = y1, group = label), size = 0.5, alpha = 0.8) + 
  geom_point( aes(x = t, y = y1, group = label), colour = "blue", alpha = 0.25) + 
  geom_line(data = DD_mu, aes(x = tgrid, y = mu1), colour = "red", size = 0.8) + 
  theme_bw()

# y2
ggplot(data = DD0) + 
  geom_line( aes(x = t, y = y2, group = label), size = 0.5, alpha = 0.8) + 
  geom_point( aes(x = t, y = y2, group = label), colour = "blue", alpha = 0.25) + 
  geom_line(data = DD_mu, aes(x = tgrid, y = mu2), colour = "red", size = 0.8) + 
  theme_bw()



#### plot covariance surfaces ####
Cov_hat <- mFPCA_out$C_hat

plot_ly() %>% 
  add_trace(
    x = 1:200,
    y = 1:200,
    z = Cov_hat,
    type = "surface"
  )


#### plot eigenfuncitons ####
eigenfunc <- mFPCA_out$Phi_hat
matplot(tgrid, eigenfunc[1:100, ], type = "l")
matplot(tgrid, eigenfunc[101:200, ], type = "l")



#### plot Xhat ####
XiEst <- mFPCA_out$Xi_hat
Xhat <- mFPCA_out$X_hat

matplot(tgrid, t(Xhat[,,1]), type = "l")
matplot(tgrid, t(Xhat[,,2]), type = "l")



## estimate vs true plot
i <- 1
k <- 1

label <- unique(DD0$label) 
index_i <- which( DD0$label == label[i] )
plot(tgrid, Xhat[i,,k], type = "l", col = "red")
lines(DD0[index_i,4], DD0[index_i,1+k])



i <- 1
k <- 2

label <- unique(DD0$label) 
index_i <- which( DD0$label == label[i] )
plot(tgrid, Xhat[i,,k], type = "l", col = "red")
lines(DD0[index_i,4], DD0[index_i,1+k])



i <- 15
k <- 2

label <- unique(DD0$label) 
index_i <- which( DD0$label == label[i] )
plot(tgrid, Xhat[i,,k], type = "l", col = "red")
lines(DD0[index_i,4], DD0[index_i,1+k])




