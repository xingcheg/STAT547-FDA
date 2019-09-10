library(fda)
library(ggplot2)
library(plotly)
library(reshape2)

######################## Problem 1 #######################

data(pinch)
dim(pinch)

##### (b) #####
pinch_mu <- apply(pinch, 1, mean)
pinch_sd <- apply(pinch, 1, sd)

r_pinch <- melt(pinch)
names(r_pinch) <- c("time", "rep", "value")
r_pinch$time <- rep(pinchtime, 20)
r_pinch$rep <- as.factor(r_pinch$rep)

## raw
ggplot(data = r_pinch) + 
  geom_line(aes(x = time, y = value, group = rep), 
              alpha = 0.5) +
  geom_line(data = data.frame(x = pinchtime, y = pinch_mu), 
            aes(x=x, y=y), colour = "red", size = 1.5) + 
  geom_line(data = data.frame(x = pinchtime, y = pinch_sd), 
            aes(x=x, y=y), colour = "blue", size = 1) +
  theme_bw()

## smoothed
ggplot(data = r_pinch) + 
  geom_smooth(aes(x = time, y = value, group = rep), 
              alpha = 0.5, span = 0.1, se = FALSE, col = "gray") +
  geom_line(data = data.frame(x = pinchtime, y = pinch_mu), 
            aes(x=x, y=y), colour = "red", size = 1.5) + 
  geom_line(data = data.frame(x = pinchtime, y = pinch_sd), 
            aes(x=x, y=y), colour = "blue", size = 1) +
  theme_bw()


##### (c) #####
pinch_Cov <- cov(t(pinch))
persp(x=pinchtime, y=pinchtime, z=pinch_Cov)
contour(x=pinchtime, y=pinchtime, z=pinch_Cov)

#plot_ly(x=pinchtime, y=pinchtime, z=pinch_Cov) %>% add_surface(
#  contours = list(
#    z = list(
#      show=TRUE,
#      usecolormap=TRUE,
#      highlightcolor="#ff0000",
#      project=list(z=TRUE)
#    )
#  )
#)



######################## Problem 2 #######################
D2 <- read.csv("/Users/apple/Desktop/ISU 2019 fall/STAT547/data/DataSets/Dow_companies_data.csv")

##### (a) #####
P1 <- D2$XOM[1]
date <- as.Date(as.character(D2$Date), format = "%m/%d/%y")

D_XOM <- data.frame(date = date, v1 = D2$XOM, 
                    v2 = 100*((D2$XOM/P1)-1))

ggplot(data = D_XOM) +
  geom_line(aes(x = date, y = v1)) + 
  ylab("stock value") + 
  ggtitle("Exxon–Mobil (XOM)") + theme_bw()

ggplot(data = D_XOM) +
  geom_line(aes(x = date, y = v2)) + 
  ylab("cumulative return function (%)") + 
  ggtitle("Exxon–Mobil (XOM)") + theme_bw()


D_XOM$v2[length(D_XOM$v2)]



##### (b) #####
D2_1 <- D2[,-1]

D2_2 <- apply(D2_1, 2, FUN = function(x){
  return( 100*(x - x[1])/x[1] )
})

r_D2_2 <- melt(D2_2)
names(r_D2_2) <- c("date", "stock", "CRF")
r_D2_2$date <- rep(date, 30)

stock_mean <- apply(D2_2, 1, mean)
stock_med <- apply(D2_2, 1, median)

DD <- data.frame(date = rep(date,2), 
                 CRF = c(stock_mean, stock_med),
                 label = rep(c("mean", "median"), each = 252))


ggplot(data = r_D2_2) + 
  geom_line(aes(x = date, y = CRF, group = stock), colour = "gray") + 
  geom_line(data = DD, aes(x = date, y = CRF, colour = label), size = 1) + 
  theme_bw()



##### (c) #####
fbplot(fit = D2_2, ylim = c(-20,90), xlab="Trading day",
       ylab = "Cumulative return")


try( fbplot(fit = D2_2, ylim = c(-20,90), prob = c(0.9, 0.6, 0.3),
       color = c(8,4,2),
       xlab="Trading day", 
       ylab = "Cumulative return"), silent = TRUE)










