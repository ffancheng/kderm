# This script is for the experiment of Spherical cap to check the H matrix from Learn Metric algorithm
rm(list = ls())
library(tidyverse)
library(scatterplot3d)
library(dimRed)
Jmisc::sourceAll(here::here("R/sources"))
set.seed(1)

# Simulate uniformly from unit sphere
n <- 10000

x <- matrix(rnorm(n * 3), n, 3)
x <- t(apply(x, 1, function(a) {a / sqrt(sum(a ^ 2))} ))
# x <-  x[x[,3] >= 0, ]; n <- nrow(x) # HALF SPHERE ONLY

theta <- acos(x[, 3]) # [0, pi], minimum angle between the position vector of given point and the +Z-axis
phi <- atan(x[, 2] / x[, 1]) # [-pi/2, pi/2], angle between the vertical half plan passing the given point and the +X-axis in anticlockwise direction, range of [0, 2pi]

phi[(x[, 2] < 0) & (x[, 1] < 0)] <- phi[(x[, 2] < 0) & (x[, 1] < 0)] + pi
phi[(x[, 2] > 0) & (x[, 1] < 0)] <- phi[(x[, 2] > 0) & (x[, 1] < 0)] - pi

y <- cbind(theta * cos(phi), theta * sin(phi))
ycap <- y[(theta < 1), ]

# Checks
cc <- rep('blue', n)
cc[theta < 1] <- 'red'
plot(y[, 1],
     y[, 2],
     col = cc,
     asp = 1,
     pch = 20)
# plot(ycap[, 1],
#      ycap[, 2],
#      col = cc[theta < 1],
#      asp = 1,
#      pch = 20)
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
caparea <- 2 * pi * 1^2 * (1 - cos(1)) # spherical cap area = 2*pi*R*h
areaball_local <- pi # area of ycap in the mapping
caparea / areaball_local # 2 * (1 - cos(1))

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
scatterplot3d(train[, 1],
              train[, 2],
              train[, 3],
              color = cc, # [theta < pi/2],
              asp = 1,
              pch = 20)
N <- nrow(train)
s <- 2
k <- 10
method <- "annIsomap" # not used if embedding y is given
annmethod <- "kdtree"
distance <- "euclidean"
treetype <- "kd"
searchtype <- "radius" # change searchtype for radius search based on `radius`, or KNN search based on `k`
radius <- 2 # radius search
# gridsize <- 20
# noutliers <- 20

suppressWarnings(metric_isomap <- metricML(x = train, s = s, 
                                           k = k, radius = radius,  # for NN search, k-NN or eps-NN depending on `searchtype`
                                           method = method, fn = fn, # for manifold learning
                                           invert.h = TRUE, 
                                           eps = 0, # for RANN::nn2(), (1+eps)*TrueDist2NN
                                           annmethod = annmethod, distance = distance, treetype = treetype,
                                           searchtype = searchtype
))
summary(metric_isomap)
fn <- metric_isomap$embedding
rmetric <- metric_isomap$rmetric # Estimated Riemannian metric
adj_matrix <- metric_isomap$adj_matrix # pairwise distance matrix within `radius`
summary(adj_matrix %>% as.vector)
# nn.idx <- metric_isomap$nn2res$nn.idx # NN index of size N

## --------------------------------------------------------
# 2) The average mean(sqrt(det(H_i))) in the spherical cap should be equal to the volume in the tangent space
# The area of the spherical cap / area of the local ball = mean squared root of the determinant of Hs (the volume of the spherical cap)
dim(rmetric)
apply(rmetric[,,theta < 1], 3, det) %>% 
  sqrt() %>% 
  # summary()
  mean()

# Tried different `n`(Line9) or `radius`(Line88) to construct the eps-NN graph to get different Laplacian matrix and the Hs are slightly different
# Use only half sphere in Line 13 to ensure that the embedding y does not include points from the "South" cap?

## Finding: The ratio of the red spherical cap area is approximately equal to the mean squared root of the determinant of the estimated pointwise H matrix from the Learn Metric algorithm
