rm(list= ls())
library(tidyverse)
library(dimRed)
library(ks)
Jmisc::sourceAll(here::here("R/sources"))

set.seed(1)
N <- 1000
d <- 1
distribution <- c("uniform", "normal", "beta")[3] # CHANGE DISTRIBUTION OF THETA
switch(distribution,
       "uniform" = {
         # 1. theta ~ U(0, pi/2)
         theta <- runif(N, 0, pi/2)
         # density
         dentheta <- dunif(theta, 0, pi/2) # density = 1/(MAX-MIN) = 0.6366198
         # plot(theta, dentheta)
       },
       "normal" = {
         # 2. theta ~ N(0, pi/20)
         mu <- 0
         sd <- pi / 20
         theta <- rnorm(N, mu, sd)
         # theta <- exp(theta)
         # hist(theta)
         dentheta <- dnorm(theta, mu, sd)
         # plot(theta, dentheta)
       },
       "beta" = {
         # 3. theta ~ Beta(a, b) \in [0,1]
         shape1 <- 2
         shape2 <- 5
         theta <- rbeta(N, shape1, shape2)
         dentheta <- dbeta(theta, shape1, shape2)
         theta <- theta * (pi / 2)
         dentheta <- dentheta / (pi / 2)
         # plot(theta, dentheta)
       })
plot(theta, dentheta)

# # NOT USED If x ~ U(0,1), theta = acos(x)
# par(mfrow=c(1,2))
# plot(x, sqrt(1-x^2), main = "Cartesian coordinates", ylim = c(0,1))
# points(x, rep(0, N), col = grey(0.5), type = "p")
# plot(x, denx, main = "U(0,1) CDF")
# # polar coordinates
# # theta <- acos(x)
# dentheta <- 1 - cos(theta)
# plot(cos(theta), sin(theta), main = "Polar coordinates", ylim = c(0,1))
# points(theta, rep(0, N), col = grey(0.5), type = "p")
# plot(theta, dentheta, main = "Polar CDF")

# Transformation
x <- cos(theta) # |dtheta/dx| = |-1/sqrt(1-x^2)|
# True density of x, f(x) = f(theta) * |dtheta/dx|
denx <- dentheta * 1 / sqrt(1 - x ^ 2)
y <- sin(theta)
plot(x, y, asp = 1)
# range(x) # [0, 1]
plot(x, denx)
# range(denx)

# KDE with fixed bandwidth
h <- ks::hpi(x, binned = TRUE) # 0.034
# h <- 0.05
fxkde <- kde(x, eval.points = x)$estimate
plot(x, fxkde)

# KDE with Riemannian metric
d <- 1
fx <- rep(0, N)
# eval.points <- x
for(i in 1:N){
  fi <- 1 / sqrt(1 - x[i] ^ 2) * dnorm(x = acos(x), mean = acos(x[i]), sd = h)
  fx <- fx + fi
}
fx <- fx / N

# Plotting
par(mfrow=c(2,2))
plot(theta, dentheta, main = "True density of theta")
plot(x, denx, main = "True density of x")
plot(x, fxkde, main = paste("Estimated f(x) with KDE, h = ", round(h, 3)))
plot(x, fx, main = "Estimated f(x) with riemannian metric")
# points(x, dentheta, col = "red", lty = "dashed")

# plot(x[order(x)][1:N*.6], fx[order(x)][1:N*.6], type = "b", cex = .5, main = "Estimated f(x)[1:.6*N] with riemannian metric")
# abline(h = dentheta, col = "red", lty = "dashed")
# text(x = 0.7, y = 0.55, paste("f(theta) = ", round(dentheta[1], 3)), col = "red")

summary(denx)
summary(fxkde)
summary(fx)
# summary(dentheta)

# rank correlation
# our estimator is better than kde
cor(denx, fx, method = "s")
cor(denx, fxkde, method = "s")
mean((denx - fx) ^ 2)
mean((denx - fxkde) ^ 2)
# mean((rank(denx) - rank(fx)) ^ 2)
# mean((rank(denx) - rank(fxkde)) ^ 2)

par(mfrow=c(2,2))
plot(rank(denx), rank(fx), main = paste("Rank correlation:", round(cor(rank(denx), rank(fx), method = "s"), 3)))
plot(rank(denx), rank(fxkde), main = paste("Rank correlation:", round(cor(rank(denx), rank(fxkde), method = "s"), 3)))
plot(denx, fx); abline(0, 1, lty = "dashed")
plot(denx, fxkde); abline(0, 1, lty = "dashed")





# Now reduce the dimension and estimate the density of the embedding
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
# E1 <- fn[,1] # rename as Ed to match the aesthetics in plot_ellipse()
xr <- diff(range(fn, na.rm = TRUE))
xextend <- 0.15
xr <- c(min(fn) - xr * xextend, max(fn) + xr * xextend)
h <- hpi(fn, binned = TRUE) # 0.0767
denks <- ks::kde(x = fn, h = h, xmin = xr[1], xmax = xr[2], gridsize = gridsize, eval.points = fn)
str(denks)
ffixed <- denks$estimate
summary(ffixed)
summary(denx)
summary(dentheta)
cor(dentheta, ffixed, method = "s") # -0.380
plot(fn, ffixed)
plot(theta, dentheta)
# plot(x, denx)

rmetric <- metric_isomap$rmetric
opt.method <- "SCALED"
riem.scale <- 1
# fisomap <- vkde(x = fn, h = Rn, gridsize = gridsize, eval.points = fn, opt.method = opt.method, riem.scale = riem.scale)
adj_matrix <- metric_isomap$adj_matrix
dim(adj_matrix)
summary(metric_isomap)
nn.idx <- metric_isomap$nn2res$nn.idx

# geodesic distance estimation
# Option 1: true geodesic distance d_g
# Option 2: approximate geodesic distance with inner product adjusted by Riem metric (x-x_i)'h_i^(-1)(x-xi)
# Option 3: approximate geodesic distance with shortest path graph distance d_G
d <- 1
fhat <- rep(0, N)
h <- h # N=100, h=0.147
# h <- 0.35
hidet <- apply(rmetric, 3, det) %>% mean()
# hidetmean <- mean(hidet)
for(i in 1:N){
  hi <- rmetric[,,i]
  # fi <- (1 - fn[i]^2) ^ (-1/2) * dnorm(x = acos(fn), mean = acos(fn[i]), sd = h) # true geodesic distance, true volume density function \theta = \sqrt{1-x^2}
  # fi <- hi ^ (-1) * dnorm(x = acos(fn), mean = acos(fn[i]), sd = h) # true geodesic distance, estimated Riem metric
  # fi <- hi ^ (-1/2) * h ^ (-1) * (2 * pi) ^ (-1/2) * exp(- adj_matrix[,i]^2 / (2 * h ^ 2)) # approximate geodesic distance with shortest path graph distance, estimated Riem metric
  # fi <- dnorm(x = fn, mean = fn[i], sd = h * sqrt(hi)) # approximate geodesic distance with inner product adjusted by Riem metric
  
  # h <- sweep(hi, 1, hidetmean / hi, "*")
  hineighbor <- mean(rmetric[,,nn.idx[i, (2:(k+1))]])
  # print(c(hi, hineighbor))
  hi <- hi * hineighbor / hi
  # hi <- hidet # gives the lowest MSE
  fi <- dnorm(x = fn, mean = fn[i], sd = h * sqrt(hi)) # approximate geodesic distance, average hi using the hi's of the neighborhood of fn[i], i.e. hi/sum(hi_neighbors)
  fhat <- fhat + fi
}
fhat <- fhat / N
plot(fn, fhat, main = paste("Estimated density of embedding", "h=", round(h,3)))
points(theta, dentheta, col = "red", lty = "dashed")
# text(x = 0.5, y = 0.5, paste("f(theta) = ", round(dentheta[1], 3)), col = "red")
summary(fhat)
cor(dentheta, fhat, method = "s")

plot(fn, ffixed)
points(theta, dentheta, col = "red", lty = "dashed")
# text(x = 0.5, y = 0.5, paste("f(theta) = ", round(dentheta[1], 3)), col = "red")
summary(ffixed)

mean((fhat - dentheta)^2)
mean((ffixed - dentheta)^2)

# How to compare rank of density estimates???
cor(dentheta, fhat, method = "s")
cor(dentheta, ffixed, method = "s")
plot(rank(dentheta), rank(fhat), main = paste("Rank correlation:", round(cor(rank(denx), rank(fhat), method = "s"), 3)))
plot(rank(dentheta), rank(ffixed), main = paste("Rank correlation:", round(cor(rank(denx), rank(ffixed), method = "s"), 3)))


## reduce N to test， N= 100， 
# true, 0.017 > 0.016
# 0.022 > 0.016
## h is optimized for kde fixed
## optimize h using AMSE for our estimator








# # # evaluate density on grid points
# # d <- 1
# # z <- NULL
# # rmetric <- metric_isomap$rmetric
# # for (i in 1:N) {
# #   hi <- as.matrix(rmetric[,,k])
# #   z <-  abind::abind(z, array(h ^ (- d / 2) * mvtnorm::dmvnorm(x = fn, mean = fn[k,], sigma = h * hi), dim = rep(N, d)), along = d + 1)  # stack array of dimension (gridsize*gridsize) with abind
# # }
# # # fhat <- rowMeans(z, dims = 2, na.rm = TRUE)
# # fhat <- rowMeans(z, na.rm = TRUE)
# # plot(fn, fhat)
# 
# 
# 
# 
# # E2 <- fn[,2]
# # prob <- c(1, 50, 95, 99) # c(1, 10, 50, 95, 99) #
# # p_hdr_isomap <- hdrscatterplot_new(E1, kde.package = "ks", levels = prob, noutliers = noutliers, label = NULL)
# # p_hdr_isomap_p <- p_hdr_isomap$p +
# #   plot_ellipse(metric_isomap, add = T, ell.no = 50, ell.size = 0,
# #                color = blues9[5], fill = blues9[1], alpha = 0.2)
# # # p_hdr_isomap
# # h_hdr_isomap <- p_hdr_isomap$den$den$h
# 
# ## ----outliers-------------------------------------------------------------------
# Rn <- metric_isomap$rmetric # array
# # fisomap <- vkde2d(x = fn[,1], y = fn[,2], h = Rn*riem.scale, gridsize = gridsize) # $x $y $z
# # fxy_isomap <- hdrcde:::interp.2d(fisomap$x, fisomap$y, fisomap$z, x0 = E1, y0 = E2)
# # plot_contour(metric_isomap, gridsize = gridsize, riem.scale = riem.scale) # estimate grid densities with vkde()
# 
# opt.method <- c("AMISE", "MEAN", "SCALED")[2]
# riem.scale <- .1
# fisomap <- vkde(x = fn, h = Rn, gridsize = gridsize, eval.points = fn, opt.method = opt.method, riem.scale = riem.scale) # 0.923
# # fisomap <- vkde(x = fn, h = Rn*riem.scale, gridsize = gridsize, eval.points = fn) # 0.64
# # # if scaling the bandwidth
# # fisomap_MEAN <- vkde(x = fn, h = Rn, gridsize = gridsize, eval.points = fn, opt.method = "MEAN") # 0.964
# # fisomap_AMISE <- vkde(x = fn, h = Rn, gridsize = gridsize, eval.points = fn, opt.method = "AMISE") # 0.941
# # h_isomap <- fisomap$H
# # str(fisomap)
# # summary(fisomap$estimate)
# # all.equal(fisomap$estimate, p_isomap$densities)
# 
# # check if vkde with grid estimate is the same as hdrcde::interp.2d
# # fixgrid_isomap <- vkde(x = fn, h = NULL, gridsize = gridsize)
# # summary(fixgrid_isomap)
# # interpden_fix <- hdrcde:::interp.2d(fixgrid_isomap$eval.points[[1]], fixgrid_isomap$eval.points[[2]], fixgrid_isomap$estimate, x0 = E1, y0 = E2)
# # all.equal(interpden_fix, fisomap$estimate)
# 
# ## ----hdroutliers----------------------------------------------------------------
# p_isomap <- plot_outlier(x = metric_isomap, gridsize = gridsize, prob = prob, riem.scale = riem.scale, f = fisomap, ell.size = 0)
# # all.equal(fxy_isomap, p_isomap$densities)
# 
# ## ----compoutlier, eval = FALSE--------------------------------------------------
# (p_isomap$p + p_hdr_isomap$p) + coord_fixed() + 
#   plot_annotation(title = "Left: variable bandwidth; Right: fixed bandwidth", theme = theme(plot.title = element_text(hjust = 0.5)))
# 
# cormethod <- c("pearson", "kendall", "spearman")[3]
# cor(preswissroll$den, p_isomap$densities, method = cormethod)
# cor(preswissroll$den, p_hdr_isomap$densities, method = cormethod)
