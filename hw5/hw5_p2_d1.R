library(dplyr) 
library(ggplot2) 
library(locfit) 
library(plotly)
library(fda)
library(reshape2)
source("/Users/apple/Desktop/ISU 2019 fall/STAT547/hw/hw5/mFPCA.R")

#### read data ####
data(gait)
Hip <- gait[, , 1]
Knee <- gait[, , 2]
time <- as.numeric( rownames(Hip) )
r_Hip <- melt(Hip)
r_Knee <- melt(Knee)
DD0 <- data.frame(label = r_Hip$Var2, y1 = r_Hip$value, y2 = r_Knee$value, 
                  t = r_Hip$Var1)
DD0 <- DD0 %>% 
  mutate(y1 = (y1 - mean(y1))/sd(y1), 
         y2 = (y2 - mean(y2))/sd(y2),
         t = (t - min(t))/diff(range(t)))



#### run mFPCA ####
mFPCA_out <- mFPCA(DD0$label, DD0$t, DD0[,c(2,3)], ngrid = 100, 
                   h_Mu = 0.1, h_Cov = 0.3, h_W = 0.5, L = 8)


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
matplot(tgrid, eigenfunc[1:100, 1:4], type = "l")
matplot(tgrid, eigenfunc[101:200, 1:4], type = "l")


#### plot Xi ####
XiEst <- mFPCA_out$Xi_hat
plot(XiEst[,1], XiEst[,2])


#### plot Xhat ####
Xhat <- mFPCA_out$X_hat

matplot(tgrid, t(Xhat[,,1]), type = "l")
matplot(tgrid, t(Xhat[,,2]), type = "l")



## estimate vs true plot
i <- 12
k <- 2
label <- unique(DD0$label) 
index_i <- which( DD0$label == label[i] )
plot(tgrid, Xhat[i,,k], type = "l", col = "red")
lines(DD0[index_i,4], DD0[index_i,1+k])



