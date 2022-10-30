## -------------------------------------------------------------------------------
## This script is for generating the 4-D manifold embedding in 100-D space and calculating its true density
## -------------------------------------------------------------------------------
# 4-D manifold: X
# 5-D semi-sphere: semisphere
# 100-D data: train, N=10000, p=100
# Gaussian mixture density: den
# True manifold density: trueden
rm(list = ls())
library(tidyverse)
library(dimRed)
library(hdrcde)
library(ks)
library(patchwork)
library(plotly)
library(mvtnorm)
Jmisc::sourceAll(here::here("R/sources"))
set.seed(1234)

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
    p2 * dmvnorm(x, rep(mu2, 4), diag(s2, 4))
}
den <- apply(X, 1, gmmdensity)

# semi-sphere radius
radius <- ceiling(max(sqrt(x1^2 + x2^2 + x3^2 + x4^2)))
radius # radius = 8
range(X)

semisphere <- cbind(X,
               x5 = sqrt(radius^2 - (x1^2 + x2^2 + x3^2 + x4^2)),
               matrix(0, N, p - 5)
)
# head(semisphere)
# semisphere %>% as_tibble() %>% summarise(x1^2 + x2^2 + x3^2 + x4^2 + x5^2)

# # distance to (0,0,0,0,0) is the same, = radius = 8
# # calculate distance to (0,0,0,0) and color them
label <- as.factor(c(rep(1, 0.99*N), rep(2, 0.01*N)))
dist2center <- sqrt(radius^2 - semisphere[,5]^2)
## Run once
# ## Plot 5-d semi-sphere
# library(tourr)
# colors <- colourvalues::colour_values(- dist2center, palette = "magma") # generate colors for locations
# pchs <- c(rep(16, 0.99*N), rep(17, 0.01*N)) # point shapes for kernels
# animate_xy(semisphere[,1:5], col = colors, pch = pchs, cex = 0.8#,
#            # axes = "bottomleft", fps = 15
#            )
# # "figures/tourr_5d_semisphere.png"
# ## Render frames as png files in a folder and convert to gif with ImageMagick
# gifpath <- "paper/figures/tourr_render/"
# tourr::render(semisphere[,1:5], grand_tour(), display_xy(col = colors, pch = pchs, cex = 0.8), "png",
#               frames = 200,
#               file.path(gifpath, "tourr_5d_semisphere-%20d.png")
#               )
## Save gif file as online supplementary file
# system("convert -delay 10 paper/figures/tourr_render/*.png paper/figures/tourr_5d_animation.gif")
## Not run
# library(gifski) # installation failed
# gif_file <- file.path('paper/figures/tourr_5d_semisphere.gif')
# save_gif(
#   makeplot(),
#          # animate_xy(semisphere[,1:5], col = colors, pch = pchs, cex = 0.8),
#          gif_file, 1280, 720, res = 144)
# utils::browseURL(gif_file)

# QR decomposition
y <- matrix(runif(p*p, 0, 1), p, p)
QR <- qr(y)
QR$rank
Q <- qr.Q(QR); #Q
train <- t(Q %*% t(semisphere))
dim(train)


## -------------------------------------------------------------------------------
# TRUE density of manifold using DC-KDE, KNN priority search
## -------------------------------------------------------------------------------
# Parameters fixed
x <- train
y <- X # given known embedding with known density
# y <- NULL
s <- 2
k <- 100 # N / 20
method <- "Isomap" # not used because y is given
annmethod <- "kdtree"
distance <- "euclidean"
treetype <- "kd"
searchtype <- "priority" # change searchtype for radius search based on `radius`, or KNN search based on `k`
radius <- 8 # the bandwidth parameter, \sqrt(\elsilon), as in algorithm. Note that the radius need to be changed for different datasets, not to increase k

gridsize <- 20
noutliers <- 20

metric_meta <- metricML(x, fn = y, s = s, k = k, radius = radius, method = method, invert.h = TRUE, eps = 0,
                        annmethod = annmethod, distance = distance, treetype = treetype,
                        searchtype = searchtype
)
# all.equal(y, metric_meta$embedding, check.attributes = F)
Rn <- metric_meta$rmetric

# Transformed from sr$den
den_4dmanifold <- den * (apply(Rn, 3, det)) ^ (.5)
summary(den_4dmanifold)
trueden <- den_4dmanifold

save(X, semisphere, train, N, p, den, trueden, k, file = paste0("data/simdata_100d_4dmanifold_N10000_trueden_k", k, ".rda"))
