# Run on HPC for N=10000, Isomap
library(tidyverse)
library(dimRed)
# library(reticulate)
library(here)
library(viridis)
library(hdrcde)
library(igraph)
# library(matrixcalc)
# library(akima)
# library(car)
library(ggforce)
library(ks)
library(patchwork)
library(plotly)
library(mvtnorm)
Jmisc::sourceAll(here::here("R/sources"))
set.seed(1234)
# scen <- as.numeric(commandArgs()[[6]])

paste("Start at:", Sys.time())

N <- 10000
p <- 100
mu1 = mu2 = 0
s1 <- 1#/(15^2)
s2 <- 2#/(15^2)
p1 <- 0.99
p2 <- 0.01
x1 <- c(rnorm(0.99*N, mu1, s1), rnorm(0.01*N, mu2, s2))
x2 <- c(rnorm(0.99*N, mu1, s1), rnorm(0.01*N, mu2, s2))
x3 <- c(rnorm(0.99*N, mu1, s1), rnorm(0.01*N, mu2, s2))
x4 <- c(rnorm(0.99*N, mu1, s1), rnorm(0.01*N, mu2, s2))
X <- cbind(x1, x2, x3, x4)

# GMM density as true meta data density
gmmdensity <- function(x) {
  p1 * dmvnorm(x, rep(mu1, 4), diag(s1, 4)) + 
    p2 *  dmvnorm(x, rep(mu2, 4), diag(s2, 4))
}
den <- apply(X, 1, gmmdensity)

# semi-sphere radius
r <- ceiling(max(sqrt(x1^2 + x2^2 + x3^2 + x4^2)))
r # radius = 15 for sigma=5
range(X)

X_new <- cbind(X, 
               x5 = sqrt(r^2 - (x1^2 + x2^2 + x3^2 + x4^2)),
               matrix(0, N, p - 5)
)
# head(X_new)
# X_new %>% as_tibble() %>% summarise(x1^2 + x2^2 + x3^2 + x4^2 + x5^2)

# distance to (0,0,0,0,0) is the same, = radius = 7
# calculate distance to (0,0,0,0) and color them
label <- as.factor(c(rep(1, 0.99*N), rep(2, 0.01*N)))
dist2center <- sqrt(r^2 - X_new[,5]^2)

# ## Plot 5-d semi-sphere, run once
# library(tourr)
# colors <- colourvalues::colour_values(- dist2center, palette = "magma") # generate colors for locations
# pchs <- c(rep(16, 0.99*N), rep(17, 0.01*N)) # point shapes for kernels
# animate_xy(X_new[,1:5], col = colors, pch = pchs, cex = 0.8,
#            axes = "bottomleft", fps=15
#            )
# # "figures/tourr_5d_semisphere.png"

y <- matrix(runif(p*p, 0, 1), p, p)
QR <- qr(y)
QR$rank
Q <- qr.Q(QR); #Q

train <- t(Q %*% t(X_new))
dim(train)



paste("ML started at:", Sys.time())

# ----parameters-----------------------------------------------------------------
# Parameters fixed
x <- train
N <- nrow(x)
s <- 5 # embedded in 5-D
k <- 20
annmethod <- "kdtree"
distance <- "euclidean"
treetype <- "kd"
searchtype <- "radius" # change searchtype for radius search based on `radius`, or KNN search based on `k`
radius <- 10 # the bandwidth parameter, \sqrt(\elsilon), as in algorithm. Note that the radius need to be changed for different datasets, not to increase k

riem.scale <- 20 # tune parameter
gridsize <- 20

## ----isomap-------------------------------------------------------------
x <- train
method <- "Isomap"
metric_isomap <- metricML(x, s = s, k = k, radius = radius, method = method, invert.h = TRUE, eps = 0,
                          # annmethod = annmethod, distance = distance, treetype = treetype, 
                          searchtype = searchtype
)
summary(metric_isomap)

fn <- metric_isomap$embedding
dim(fn)

## ----Fixed bandwidth density estimates with ks::kde-------------------------
fixden_isomap <- vkde(x = fn, h = NULL, eval.points = fn, binned = FALSE) # gridsize = gridsize
summary(fixden_isomap$estimate)
fixden_isomap$H

## ----VKDE-------------------------------------------------------------------
Rn <- metric_isomap$rmetric # array
tictoc::tic()
fisomap <- vkde(x = fn, h = Rn*riem.scale, eval.points = fn)
tictoc::toc()
summary(fisomap$estimate)

# ----compare with true density----------------------------------------------
png(filename=paste0("figures/compareden_5d_N", N, "_", method, "_riem", format(riem.scale, decimal.mark = "_"), ".png"), width = 10, height = 6, dpi = 300)
par(mfrow=c(1,2))
plot(den, fisomap$estimate, main = paste("Variable bandwidth correlation:", round(cor(den, fisomap$estimate), 3)))
plot(den, fixden_isomap$estimate, main = paste("Fixed bandwidth correlation:", round(cor(den, fixden_isomap$estimate), 3)))
dev.off()
# dev.print(pdf, 'filename.pdf')

cor(den, fisomap$estimate)
cor(den, fixden_isomap$estimate)

save(method, metric_isomap, fixden_isomap, fisomap, train, den, file = paste0("figures/compareden_5d_N", N, "_", method, "_riem", format(riem.scale, decimal.mark = "_"), ".rda"))

paste("End at:", Sys.time())
