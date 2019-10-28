## problem 2

library (SemiPar)
library(locfit) 
library(plotly)

data(scallop)

## plot raw data
p0 <- plot_ly() %>% 
  add_trace(
  data = scallop,
  x = ~latitude, 
  y = ~longitude, 
  z = ~tot.catch,
  type = "scatter3d",
  mode = "markers",
  size = 3
  )


## plot log data
eps <- 0.1
p <- plot_ly() %>% 
  add_trace(
  data = scallop,
  x = ~latitude, 
  y = ~longitude, 
  z = ~log( tot.catch + eps),
  type = "scatter3d",
  mode = "markers",
  size = 3
)


## locpol fit using gcv
m <- 50

## gcv
nnCandidate <- seq(0.1, 1, 0.1)
gcvScores <- sapply( 
  nnCandidate, 
  function(nn) { 
    gcv( log(tot.catch + eps) ~ lp(latitude, longitude, deg=1, nn=nn),  
         scallop, ev=lfgrid(mg=m))['gcv'] 
  } 
) 

nnGCV <- nnCandidate[which.min(gcvScores)] 

## fit
res <- locfit( log(tot.catch + eps) ~ lp(latitude, longitude, deg=1, nn = nnGCV),
               scallop, ev=lfgrid(mg=m)) 


xGrid <- seq(min(scallop$latitude), max(scallop$latitude), length.out=m) 
yGrid <- seq(min(scallop$longitude), max(scallop$longitude), length.out=m) 
z <- matrix(predict(res), m, m, byrow = TRUE)



p1 <- p %>% 
  add_trace(
  x = xGrid,
  y = yGrid,
  z = z,
  type = "surface"
)


p0
p
p1



## problem 3
library(fdapace)
yeast <- read.table(file = "/Users/apple/Desktop/ISU 2019 fall/STAT547/data/yeast.txt")
time <- as.numeric( sub("alpha", "", names(yeast)) )

## original data
X <- as.matrix(yeast)
matplot(time, t(X), type='l', col=c(rep(2,44), rep(4,45)), 
        xlab = "time", ylab = "yeast") 

data_o <- MakeFPCAInputs(tVec = time, yVec = X)
res_o <- FPCA(data_o$Ly, data_o$Lt) 
plot(res_o)


## smoothed data
n <- nrow(X)
m <- 50
X_smooth <- matrix(0, ncol = m, nrow = n)
tGrid <- seq(range(time)[1], range(time)[2], length.out = m)
for (i in 1:n){
  fit <- locfit( X[i,] ~ lp(time, deg=3), ev=lfgrid(mg=m)) 
  X_smooth[i,] <- predict(fit)
}

matplot(tGrid, t(X_smooth), type='l', col=c(rep(2,44), rep(4,45)), 
        xlab = "time", ylab = "yeast") 

data_s <- MakeFPCAInputs(tVec = tGrid, yVec = X_smooth)
res_s <- FPCA(data_s$Ly, data_s$Lt) 
plot(res_s)


## derivative
X_deriv <- matrix(0, ncol = m, nrow = n)
for (i in 1:n){
  fit <- locfit( X_smooth[i,] ~ lp(tGrid, deg=4, nn=0.2), ev=lfgrid(mg=m), deriv=1) 
  X_deriv[i,] <- predict(fit)
}

matplot(tGrid, t(X_deriv), type='l', col=c(rep(2,44), rep(4,45)), 
        xlab = "time", ylab = "yeast") 

data_d <- MakeFPCAInputs(tVec = tGrid, yVec = X_deriv)
res_d <- FPCA(data_d$Ly, data_d$Lt) 
plot(res_d)




