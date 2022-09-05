# This script is for the experiment of Spherical cap to check the H matrix from Learn Metric algorithm
rm(list = ls())
library(tidyverse)
library(scatterplot3d)
library(dimRed)
library(kableExtra)
Jmisc::sourceAll(here::here("R/sources"))
set.seed(1)

# Simulate uniformly from unit sphere
n <- 2000

x <- matrix(rnorm(n * 3), n, 3)
x <- t(apply(x, 1, function(a) {a / sqrt(sum(a ^ 2))} ))
x <-  x[x[,3] >= 0, ]; n <- nrow(x) # HALF SPHERE ONLY

theta <- acos(x[, 3]) # [0, pi], minimum angle between the position vector of given point and the +Z-axis
phi <- atan(x[, 2] / x[, 1]) # atan() returns [-pi/2, pi/2]; phi is the angle between the vertical half plane passing the given point and the +X-axis in anticlockwise direction, range of [0, 2pi]
# summary(phi)
## Correct phi
# if {x>0} v=atan(y/x);
# if {y>=0 & x<0} v=pi+atan(y/x);
# if {y<0 & x<0} v=-pi+atan(y/x);
phi[(x[, 2] < 0) & (x[, 1] < 0)] <- phi[(x[, 2] < 0) & (x[, 1] < 0)] - pi
phi[(x[, 2] > 0) & (x[, 1] < 0)] <- phi[(x[, 2] > 0) & (x[, 1] < 0)] + pi

theta0 <- .5
y <- cbind(theta * cos(phi), theta * sin(phi))
ycap <- y[(theta < theta0), ]

## To compare with python megaman implementation
# library(feather)
# df <- data.frame(x); colnames(df) = c("x", "y", "z")
# write_feather(df, "python/sphericalcap_x.feather")
# write_feather(data.frame(y), "python/sphericalcap_y.feather")
# write_feather(data.frame(ycap), "python/sphericalcap_ycap.feather")

# Checks
cc <- rep('blue', n)
cc[theta < theta0] <- 'black'
plot(y[, 1],
     y[, 2],
     col = cc,
     asp = 1,
     pch = 20)
plot(ycap[, 1],
     ycap[, 2],
     col = cc[theta < theta0],
     asp = 1,
     pch = 20)
scatterplot3d(x[, 1],
              x[, 2],
              x[, 3],
              color = cc,
              asp = 1,
              pch = 20)

# sx <- cbind(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)) # r=1
# scatterplot3d(sx[, 1],
#               sx[, 2],
#               sx[, 3],
#               color = cc,
#               asp = 1,
#               pch = 20)


# https://mathworld.wolfram.com/SphericalCap.html
# https://www.gradplus.pro/gate-questions/what-is-spherical-coordinate-system/
# ratio of spherical cap area
caparea <- 2 * pi * 1 * 1 * (1 - cos(theta0)) # spherical cap area = 2*pi*R*h, with R=1
areaball_local <- pi * (1 * theta0) ^ 2 # area of ycap in the tangent space, square with radius 1* theta0, the geodesic distance
(arearatio <- caparea / areaball_local) # 2 * (1 - cos(theta0)) / (theta0^2)

# # spherical cap volume
# capvol <- 1 / 3 * pi * (1 ^ 3) * (2 - 3 * cos(1) + cos(1)^3)

## --------------------------------------------------------
## Next Steps

## 1) Treat y as an output embedding and x as an input. Run the LearnMetric algorithm to get a H at each y point
## 2) Find Average root determinant of H in the spherical cap

## --------------------------------------------------------
# 1) Input: high-dim x; Output: low-dim embedding y
# Now reduce the dimension and estimate the density of the embedding
train <- x # half sphere
fn <- y # known embedding
# fn <- NULL # unknown embedding
scatterplot3d(train[, 1],
              train[, 2],
              train[, 3],
              color = cc, # [theta < pi/2],
              asp = 1,
              pch = 20)
N <- nrow(train)
s <- 2
k <- N-1
eps <- 0
method <- "annIsomap" # not used if embedding y is given
annmethod <- "kdtree"
distance <- "euclidean"
treetype <- "kd"
searchtype <- "priority" # change searchtype for radius search based on `radius`, or "standard" / "priority" KNN search based on `k`
radius <- 2 # for radius search
# gridsize <- 20
# noutliers <- 20


# a <- rep(NA, n-2)
# for(k in 2:(n-1)) {
  suppressWarnings(metric_isomap <- metricML(x = train, s = s, 
                                             k = k, radius = radius,  # for NN search, k-NN or eps-NN depending on `searchtype`
                                             bandwidth = 0.4, # for weighted NN matrix
                                             method = method, fn = fn, # for manifold learning
                                             invert.h = TRUE, 
                                             eps = 0, # for RANN::nn2(), (1+eps)*TrueDist2NN
                                             annmethod = annmethod, distance = distance, treetype = treetype,
                                             searchtype = searchtype
  ))
  summary(metric_isomap)
  fn <- metric_isomap$embedding
  rmetric <- metric_isomap$rmetric # Estimated Riemannian metric
  weight_matrix <- metric_isomap$weighted_matrix # pairwise distance matrix within `radius`
  summary(metric_isomap$adj_matrix %>% as.vector)
  # nn.idx <- metric_isomap$nn2res$nn.idx # NN index of size N
  
  ## --------------------------------------------------------
  # 2) The average mean(sqrt(det(H_i))) in the spherical cap should be equal to the volume in the tangent space
  # The area of the spherical cap / area of the local ball = mean squared root of the determinant of Hs (the volume of the spherical cap)
  # dim(rmetric)
  # a[k-1] <- 
    # apply(rmetric[,,theta < theta0], 3, det) %>% 
    # sqrt() %>% 
    # mean()
# }
# 
# plot(a)

a <- apply(rmetric[,,theta < theta0], 3, det) %>% sqrt()
paste("Estimated ratio of area tangent space/spherical cap:", mean(a)); paste("True ratio of area tangent space/spherical cap:", areaball_local / caparea)
paste("Estimated ratio of area spherical cap/tangent space:", 1 / mean(a)); paste("True ratio of area spherical cap/tangent space:", caparea / areaball_local)

# Tried different `n`(Line 10) or `theta0`(Line 26) or `radius`(Line 103) to construct the eps-NN graph to get different Laplacian matrix and the Hs are slightly different
# Use only half sphere in Line 13 to ensure that the embedding y does not include points from the "South" cap?

## Finding: The ratio of the red spherical cap area is approximately equal to the mean squared root of the determinant of the estimated pointwise H matrix from the Learn Metric algorithm




## --------------------------------------------------------
## Verify volume density function form
## --------------------------------------------------------
# Now pick one point, the north point p=(0,0,1)
# Use ML to learn the embedding, and locate the embedded north point with new coordinates f(p) adn H(p). Note that any ML embedding could be used.
# To transform coordinates, the volume density function is used, i.e. f\prime(q)=det(H(p))^{-1/2}f(q) for all q.
# TODO: verify the right form for the volume density function, \sqrt(|H(p)|) and \sqrt(|\bar{Hi}|)
train <- rbind(x, c(0,0,1))
fn <- NULL # rbind(y, c(0,0)) # NULL if no pre-computed embedding
N <- nrow(train)
s <- 2
k <- N-1
eps <- 0
method <- "annIsomap" # not used if embedding y is given
annmethod <- "kdtree"
distance <- "euclidean"
treetype <- "kd"
searchtype <- "priority" # change searchtype for radius search based on `radius`, or "standard" / "priority" KNN search based on `k`
radius <- 2 # radius search
# gridsize <- 20
# noutliers <- 20

suppressWarnings(metric_isomap <- metricML(x = train, s = s, 
                                           k = k, radius = radius,  # for NN search, k-NN or eps-NN depending on `searchtype`
                                           bandwidth = 0.4, # for weighted NN matrix
                                           method = method, fn = fn, # for manifold learning
                                           invert.h = TRUE, 
                                           eps = 0, # for RANN::nn2(), (1+eps)*TrueDist2NN
                                           annmethod = annmethod, distance = distance, treetype = treetype,
                                           searchtype = searchtype
))
summary(metric_isomap)
fn <- metric_isomap$embedding
rmetric <- metric_isomap$rmetric
adj_matrix <- metric_isomap$adj_matrix
nn.idx <- metric_isomap$nn2res$nn.idx

# Embedded north point p=(0,0,1)
# Coordinates in the tangent space TpM is (0,0)
fn[N,]
rmetric[,,N]
# adj_matrix[N,]

# Metric embedding of p: locally isometric embedding
fn_metric <- det(rmetric[,,N]) ^ (-1/2) * fn
fn_metric[N,]

plot(fn_metric, asp = 1, col = c(cc, "red"), xlim = c(-1,1), ylim=c(-1,1))
plot(rbind(y,c(0,0)), asp = 1, col = c(cc, "red"), xlim = c(-1,1), ylim=c(-1,1))
mean(sum((rbind(y, c(0,0)) - fn_metric)^2))

# Other forms
det(rmetric[,,N]) ^ (-1/2) * fn[N,]
sqrt(det(rmetric[,,N]) * mean(apply(rmetric, 3, det)[nn.idx[N,]]) ) * fn[N,]
sqrt(det(rmetric[,,N]) / mean(apply(rmetric, 3, det)[nn.idx[N,]]) ) * fn[N,]
sqrt(mean(apply(rmetric, 3, det)[nn.idx[N,]]) / det(rmetric[,,N]) ) * fn[N,]



## TODO: check local coordinate change 
# g_q = h_p ^ {-.5} %*% h_q %*% h_p ^ {-.5}
rmetric_det <- apply(rmetric, 3, det)
r <- 1
fn_dc <- matrix(NA, N, 2)
theta <- rep(NA, N)
for(i in 1:N){
  bindex <- which(adj_matrix[i,] <= r)
  # theta[i] <- sqrt(rmetric_det[i])
  # theta[i] <- sqrt(rmetric_det[i]) ^ (-1)
  # theta[i] <- sqrt(mean(rmetric_det[bindex]) * rmetric_det[i] ) ^ (-1)
  # theta[i] <- sqrt(mean(rmetric_det[bindex]) / rmetric_det[i] ) ^ (-1)
  theta[i] <- sqrt( mean(rmetric_det[bindex]) / rmetric_det[i] ) # used in the paper, (\frac{|\det \pmb{H}(\pmb{y}_i)|}{|\det \pmb{H}(\pmb{p})|} )^{1/2}
  fn_dc[i,] <- fn[i,] * theta[i]
}
par(mfrow=c(1,3))
plot(rbind(y, c(0,0)), asp = 1, col = c(cc, "red"), xlim = c(-2, 2), ylim=c(-2, 2))
plot(fn, asp = 1, col = c(cc, "red"), xlim = c(-2, 2), ylim=c(-2, 2))
plot(fn_dc, asp = 1, col = c(cc, "red"), xlim = c(-2, 2), ylim=c(-2, 2))






## Assume that the distortion-corrected coordinates gives another embedding
# Now the problem is to compare which embedding is better
# fn or fn_dc
# compared with the logarithm map around p=c(0,0,1), rbind(y, c(0,0)

# cap_index <- c(which(cc == "black"), N) ## to check only the cap area
cap_index <- 1:N ## to check the half sphere

# Three objects to be compared
f_log <- rbind(y, c(0,0))[cap_index,]
f_ml <- fn[cap_index,]
f_dc <- fn_dc[cap_index,]

# summary(f_log))
# summary(f_ml)
# summary(f_dc)

# Sum of squared error
sum((f_log - f_dc) ^ 2)
sum((f_log - f_ml) ^ 2)

# procrustes measure
vegan::procrustes((f_log), f_dc)
vegan::procrustes((f_log), f_ml)

# embedding qualities
dr_quality(f_log, f_dc, K = 20)
# $quality
# M_T       M_C      LCMC      Qnx         W_n       W_nu       Rnx
# 0.9991575 0.9993159 0.8638992 0.884391 0.001838109 0.00159078 0.8819724

dr_quality(f_log, f_ml, K = 20)
# $quality
# M_T       M_C      LCMC       Qnx         W_n        W_nu       Rnx
# 0.9989322 0.9979771 0.8183004 0.8387922 0.002183664 0.002038867 0.8354197


measures <- tibble(sse = c(sum((f_log - f_dc) ^ 2),
               sum((f_log - f_ml) ^ 2)),
       proc = c(vegan::procrustes((f_log), f_dc)$ss,
                vegan::procrustes((f_log), f_ml)$ss),
       rbind(dr_quality(f_log, f_dc, K = 20)$quality,
             dr_quality(f_log, f_ml, K = 20)$quality)
       )
kable(measures, booktabs = TRUE, digits = 5, escape = FALSE) %>%
  kable_styling(latex_options = "scale_down") %>%
  kable_paper(full_width = FALSE) 
