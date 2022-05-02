rm(list= ls())
library(tidyverse)
library(dimRed)
# library(viridis)
# library(hdrcde)
library(ks)
# library(patchwork)
# # library(copula)
# library(plotly)
# library(kableExtra)
Jmisc::sourceAll(here::here("R/sources"))

set.seed(1)
N <- 10000
d <- 2
theta <- runif(N, 0, pi/2)
# theta <- rnorm(N, pi / 4, pi / 20)
# density
dentheta <- dunif(theta, 0, pi/2) # density = 1/(MAX-MIN) = 0.6366198
x <- cos(theta)
y <- sin(theta)
plot(x, y, asp = 1)
# range(x) # [0, 1]
denx <- 2 / pi * (1 / sqrt(1 - x ^ 2))
plot(x, denx)
range(denx)


h <- ks::hpi(x, binned = TRUE) # 0.034
# h <- 0.05
fxkde <- kde(x, eval.points = x)$estimate
plot(x, fxkde)

d <- 1
fx <- rep(0, N)
# eval.points <- x
for(i in 1:N){
  # fi <- sqrt(1 - x[i] ^ 2) * dmvnorm(x = eval.points, mean = x[i], sigma = diag(h, N))
  # fx <- sum(fx, fi)
  fi <- 1 / sqrt(1 - x[i] ^ 2) * dnorm(x = acos(x), mean = acos(x[i]), sd = h)
  fx <- fx + fi
}
fx <- fx / N

par(mfrow=c(2,2))
plot(x, denx)
plot(x, fx)
plot(x, fxkde)
plot(x, dentheta)

cor(denx, fx, method = "s")
cor(denx, fxkde, method = "s")

hist(x, probability = T)



# evaluate density on grid points
fhat <- NULL
for (k in 1:N) {
  # hk <- h[,,k]
  z <-  abind::abind(z, array(mvtnorm::dmvnorm(x = eval.points, mean = x[k,], sigma = hk), dim = gridsize), along = d + 1)  # stack array of dimension (gridsize*gridsize) with abind
}
z <- rowMeans(z, dims = 2, na.rm = TRUE)









# # If x ~ U(0,1), theta = acos(x)
# par(mfrow=c(1,2))
# plot(x, sqrt(1-x^2), main = "Cartesian coordinates", ylim = c(0,1))
# points(x, rep(0, N), col = grey(0.5), type = "p")
# plot(x, denx, main = "U(0,1) CDF")
# 
# # polar coordinates
# # theta <- acos(x)
# dentheta <- 1 - cos(theta)
# plot(cos(theta), sin(theta), main = "Polar coordinates", ylim = c(0,1))
# points(theta, rep(0, N), col = grey(0.5), type = "p")
# plot(theta, dentheta, main = "Polar CDF")


train <- data.frame(x = x, y = y)
# Parameters fixed
# x <- train
N <- nrow(train)
s <- 1
k <- 10
method <- "Isomap"
annmethod <- "kdtree"
distance <- "euclidean"
treetype <- "kd"
searchtype <- "radius" # change searchtype for radius search based on `radius`, or KNN search based on `k`
radius <- 8 # the bandwidth parameter, \sqrt(\elsilon), as in algorithm. Note that the radius need to be changed for different datasets, not to increase k

gridsize <- 20
noutliers <- 20
# riem.scale <- .1 # .1 # scale Riemmanian

# ISOMAP

## ---- message=FALSE-------------------------------------------------------------
metric_isomap <- metricML(x = train, s = s, k = k, radius = 2, method = method, invert.h = TRUE, eps = 0,
                          # annmethod = annmethod, distance = distance, treetype = treetype, 
                          searchtype = searchtype
)
# summary(metric_isomap)

## -------------------------------------------------------------------------------
# fixed bandwidth
fn <- metric_isomap$embedding
E1 <- fn[,1] # rename as Ed to match the aesthetics in plot_ellipse()
xr <- diff(range(fn, na.rm=TRUE))
xextend <- 0.15
xr <- c(min(x)-xr*xextend,max(x)+xr*xextend)
h <- hpi(fn, binned = TRUE) # 0.0767
denks <- ks::kde(x = fn, h = h, xmin = xr[1], xmax = xr[2], gridsize = gridsize, eval.points = fn)
str(denks)
denfixed <- denks$estimate
range(denfixed)
range(denx)
cor(denx, denfixed, method = "s") # -0.380

Rn <- metric_isomap$rmetric
opt.method <- "SCALED"
riem.scale <- 1
# fisomap <- vkde(x = fn, h = Rn, gridsize = gridsize, eval.points = fn, opt.method = opt.method, riem.scale = riem.scale)














# E2 <- fn[,2]
# prob <- c(1, 50, 95, 99) # c(1, 10, 50, 95, 99) #
# p_hdr_isomap <- hdrscatterplot_new(E1, kde.package = "ks", levels = prob, noutliers = noutliers, label = NULL)
# p_hdr_isomap_p <- p_hdr_isomap$p +
#   plot_ellipse(metric_isomap, add = T, ell.no = 50, ell.size = 0,
#                color = blues9[5], fill = blues9[1], alpha = 0.2)
# # p_hdr_isomap
# h_hdr_isomap <- p_hdr_isomap$den$den$h

## ----outliers-------------------------------------------------------------------
Rn <- metric_isomap$rmetric # array
# fisomap <- vkde2d(x = fn[,1], y = fn[,2], h = Rn*riem.scale, gridsize = gridsize) # $x $y $z
# fxy_isomap <- hdrcde:::interp.2d(fisomap$x, fisomap$y, fisomap$z, x0 = E1, y0 = E2)
# plot_contour(metric_isomap, gridsize = gridsize, riem.scale = riem.scale) # estimate grid densities with vkde()

opt.method <- c("AMISE", "MEAN", "SCALED")[2]
riem.scale <- .1
fisomap <- vkde(x = fn, h = Rn, gridsize = gridsize, eval.points = fn, opt.method = opt.method, riem.scale = riem.scale) # 0.923
# fisomap <- vkde(x = fn, h = Rn*riem.scale, gridsize = gridsize, eval.points = fn) # 0.64
# # if scaling the bandwidth
# fisomap_MEAN <- vkde(x = fn, h = Rn, gridsize = gridsize, eval.points = fn, opt.method = "MEAN") # 0.964
# fisomap_AMISE <- vkde(x = fn, h = Rn, gridsize = gridsize, eval.points = fn, opt.method = "AMISE") # 0.941
# h_isomap <- fisomap$H
# str(fisomap)
# summary(fisomap$estimate)
# all.equal(fisomap$estimate, p_isomap$densities)

# check if vkde with grid estimate is the same as hdrcde::interp.2d
# fixgrid_isomap <- vkde(x = fn, h = NULL, gridsize = gridsize)
# summary(fixgrid_isomap)
# interpden_fix <- hdrcde:::interp.2d(fixgrid_isomap$eval.points[[1]], fixgrid_isomap$eval.points[[2]], fixgrid_isomap$estimate, x0 = E1, y0 = E2)
# all.equal(interpden_fix, fisomap$estimate)

## ----hdroutliers----------------------------------------------------------------
p_isomap <- plot_outlier(x = metric_isomap, gridsize = gridsize, prob = prob, riem.scale = riem.scale, f = fisomap, ell.size = 0)
# all.equal(fxy_isomap, p_isomap$densities)

## ----compoutlier, eval = FALSE--------------------------------------------------
(p_isomap$p + p_hdr_isomap$p) + coord_fixed() + 
  plot_annotation(title = "Left: variable bandwidth; Right: fixed bandwidth", theme = theme(plot.title = element_text(hjust = 0.5)))

cormethod <- c("pearson", "kendall", "spearman")[3]
cor(preswissroll$den, p_isomap$densities, method = cormethod)
cor(preswissroll$den, p_hdr_isomap$densities, method = cormethod)
