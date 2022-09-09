# This script is for the experiment of Spherical cap to check the H matrix from Learn Metric algorithm
rm(list = ls())
library(tidyverse)
library(scatterplot3d)
library(dimRed)
library(kableExtra)
Jmisc::sourceAll(here::here("R/sources"))
set.seed(1)

# Simulate from unit sphere
n <- 3000
## --------------------------------------------------------
# x <- matrix(rnorm(n * 3), n, 3)
# x <- t(apply(x, 1, function(a) {a / sqrt(sum(a ^ 2))} ))
# x <-  x[x[,3] >= 0, ]; n <- nrow(x) # HALF SPHERE ONLY
# 
# theta <- acos(x[, 3]) # [0, pi], minimum angle between the position vector of given point and the +Z-axis
# phi <- atan(x[, 2] / x[, 1]) # atan() returns [-pi/2, pi/2]; phi is the angle between the vertical half plane passing the given point and the +X-axis in anticlockwise direction, range of [0, 2pi]
# # summary(phi)
# ## Correct phi
# # if {x>0} v=atan(y/x);
# # if {y>=0 & x<0} v=pi+atan(y/x);
# # if {y<0 & x<0} v=-pi+atan(y/x);
# phi[(x[, 2] < 0) & (x[, 1] < 0)] <- phi[(x[, 2] < 0) & (x[, 1] < 0)] - pi
# phi[(x[, 2] > 0) & (x[, 1] < 0)] <- phi[(x[, 2] > 0) & (x[, 1] < 0)] + pi
# 
# theta0 <- .5
# y <- cbind(theta * cos(phi), theta * sin(phi))
# ycap <- y[(theta < theta0), ]
# 
# ## To compare with python megaman implementation
# # library(feather)
# # df <- data.frame(x); colnames(df) = c("x", "y", "z")
# # write_feather(df, "python/sphericalcap_x.feather")
# # write_feather(data.frame(y), "python/sphericalcap_y.feather")
# # write_feather(data.frame(ycap), "python/sphericalcap_ycap.feather")
# 
# # Checks
# cc <- rep('blue', n)
# cc[theta < theta0] <- 'black'
# plot(y[, 1], y[, 2], col = cc,
#      asp = 1, pch = 20)
# plot(ycap[, 1], ycap[, 2], col = cc[theta_cap],
#      asp = 1, pch = 20)
# scatterplot3d(x[, 1], x[, 2], x[, 3], color = cc,
#               asp = 1, pch = 20)


## --------------------------------------------------------
# OR: Simulate theta and phi uniformly and transform polar coordinates
theta <- runif(n, 0, pi)
phi <- runif(n, 0, 2 * pi)
p_A_theta_phi <- 1 / (2 * pi ^ 2)
den_S2 <- p_A_theta_phi / sin(theta) ## IMPORTANT
sx <- cbind(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)) # r=1
cc <- rep('blue', n)
theta0 <- .5
cc[theta < theta0] <- 'black'
half_index <- sx[ ,3] >= 0
scatterplot3d(sx[, 1],
              sx[, 2],
              sx[, 3],
              color = cc,
              asp = 1,
              pch = 20)

x <-  sx[half_index, ]; n <- nrow(x) # HALF SPHERE ONLY
theta_cap <- (theta < theta0) & (sx[ ,3] >= 0) # index in whole sphere
y <- cbind(theta * cos(phi), theta * sin(phi))
ycap <- y[theta_cap, ]
y <- y[half_index, ]
plot(ycap[, 1], ycap[, 2], col = cc[theta_cap],
     asp = 1, pch = 20)
cap_half_index <- which(theta[half_index] < theta0) # index in half sphere
cc <- cc[half_index]
plot(y[, 1], y[, 2], col = cc,
     asp = 1, pch = 20)
scatterplot3d(x[, 1], x[, 2], x[, 3], color = cc,
              asp = 1, pch = 20)

## --------------------------------------------------------
## True density of half sphere and cap area
den_half_sphere <- den_S2[half_index]
den_cap <- den_S2[theta_cap]
# summary(den_cap)

## --------------------------------------------------------
## Ratio of spherical cap area
# https://mathworld.wolfram.com/SphericalCap.html
# https://www.gradplus.pro/gate-questions/what-is-spherical-coordinate-system/
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
    # apply(rmetric[,,cap_half_index], 3, det) %>% 
    # sqrt() %>% 
    # mean()
# }
# 
# plot(a)

a <- apply(rmetric[,,cap_half_index], 3, det) %>% sqrt()
paste("Estimated ratio of area tangent space/spherical cap:", mean(a)); paste("True ratio of area tangent space/spherical cap:", areaball_local / caparea)
paste("Estimated ratio of area spherical cap/tangent space:", 1 / mean(a)); paste("True ratio of area spherical cap/tangent space:", caparea / areaball_local)

# Tried different `n`(Line 10) or `theta0`(Line 26) or `radius`(Line 103) to construct the eps-NN graph to get different Laplacian matrix and the Hs are slightly different
# Use only half sphere in Line 13 to ensure that the embedding y does not include points from the "South" cap?

## Finding: The ratio of the red spherical cap area is approximately equal to the mean squared root of the determinant of the estimated pointwise H matrix from the Learn Metric algorithm




## --------------------------------------------------------
## Verify volume density function form by cap area
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

## --------------------------------------------------------
# Metric embedding of p: locally isometric embedding
# fn_metric = hn ^ (-1/2) %*% fn # WRONG
hp <- expm::sqrtm(solve(rmetric[,,N])) # only using the Riemannian metric for the north point (0,0,1)
fn_metric <- t(apply(fn, 1, function(x) t(hp %*% x)))
fn_metric[N,]

par(mfrow=c(1,3))
plot(rbind(y,c(0,0)), asp = 1, col = c(cc, "red"), xlim = c(-2, 2), ylim = c(-2, 2), main = "Logarithm map")
plot(fn_metric, asp = 1, col = c(cc, "red"), xlim = c(-2, 2), ylim = c(-2, 2), main = "Metric embedding")
plot(fn, asp = 1, col = c(cc, "red"), xlim = c(-2, 2), ylim = c(-2, 2), main = "ISOMAP embedding")

# mean(sum((rbind(y, c(0,0)) - fn_metric)[c(which(cc == "black"), N)]^2))
# mean(sum((rbind(y, c(0,0)) - fn)[c(which(cc == "black"), N)]^2))
# # Other forms
# det(rmetric[,,N]) ^ (-1/2) * fn[N,]
# sqrt(det(rmetric[,,N]) * mean(apply(rmetric, 3, det)[nn.idx[N,]]) ) * fn[N,]
# sqrt(det(rmetric[,,N]) / mean(apply(rmetric, 3, det)[nn.idx[N,]]) ) * fn[N,]
# sqrt(mean(apply(rmetric, 3, det)[nn.idx[N,]]) / det(rmetric[,,N]) ) * fn[N,]

## --------------------------------------------------------
# Compare the cap area
# only rmetric[,,N] is used, so isometric only locally around north point
cap_index <- c(which(cc == "black"), N)
# cap_index <- 1:N
f_log <- rbind(y, c(0,0))[cap_index,]
f_ml <- fn[cap_index,]
f_metric <- fn_metric[cap_index,]
measures_metric <- tibble(SSE = c(sum((f_log - f_metric) ^ 2),
                           sum((f_log - f_ml) ^ 2)),
                   Proc = c(vegan::procrustes(f_log, f_metric)$ss,
                            vegan::procrustes(f_log, f_ml)$ss),
                   rbind(dr_quality(f_log, f_metric, K = 20)$quality,
                         dr_quality(f_log, f_ml, K = 20)$quality)
)
kable(measures_metric, booktabs = TRUE, digits = 5, escape = FALSE) %>%
  kable_styling(latex_options = "scale_down") %>%
  kable_paper(full_width = FALSE) 
# f_metric is NOT improving
par(mfrow=c(1,1))
plot(f_log, xlim = c(-1, 1), ylim = c(-1, 1), main = "Spherical cap embeddings")
points(f_ml, col = "red")
points(f_metric, col = "green")



## --------------------------------------------------------
## Check local coordinate change 
# g_q = h_p ^ {-.5} %*% h_q %*% h_p ^ {-.5}
# |g_q| = |h_q| / |h_p|
rmetric_det <- apply(rmetric, 3, det)
r <- 1
fn_dc <- matrix(NA, N, 2)
theta <- rep(NA, N)
for(i in 1:N){
  bindex <- which(adj_matrix[i,] <= r)
  # theta[i] <- sqrt(rmetric_det[i])
  # theta[i] <- sqrt(rmetric_det[i]) ^ (-1)
  theta[i] <- sqrt(mean(rmetric_det[bindex]) * rmetric_det[i] )
  # theta[i] <- sqrt(mean(rmetric_det[bindex]) / rmetric_det[i] ) ^ (-1)
  # theta[i] <- sqrt( mean(rmetric_det[bindex]) / rmetric_det[i] ) # used in the paper, (\frac{|\det \pmb{H}(\pmb{y}_i)|}{|\det \pmb{H}(\pmb{p})|} )^{1/2}
  fn_dc[i,] <- fn[i,] * theta[i]
}
par(mfrow=c(1,3))
plot(rbind(y, c(0,0)), asp = 1, col = c(cc, "red"), xlim = c(-2, 2), ylim=c(-2, 2))
plot(fn_dc, asp = 1, col = c(cc, "red"), xlim = c(-2, 2), ylim=c(-2, 2))
plot(fn, asp = 1, col = c(cc, "red"), xlim = c(-2, 2), ylim=c(-2, 2))






## --------------------------------------------------------
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

# # summary(f_log))
# # summary(f_ml)
# # summary(f_dc)
# 
# # Sum of squared error
# sum((f_log - f_dc) ^ 2)
# sum((f_log - f_ml) ^ 2)
# 
# # procrustes measure
# vegan::procrustes((f_log), f_dc)
# vegan::procrustes((f_log), f_ml)
# 
# # embedding qualities
# dr_quality(f_log, f_dc, K = 20)
# # $quality
# # M_T       M_C      LCMC      Qnx         W_n       W_nu       Rnx
# # 0.9991575 0.9993159 0.8638992 0.884391 0.001838109 0.00159078 0.8819724
# 
# dr_quality(f_log, f_ml, K = 20)
# # $quality
# # M_T       M_C      LCMC       Qnx         W_n        W_nu       Rnx
# # 0.9989322 0.9979771 0.8183004 0.8387922 0.002183664 0.002038867 0.8354197


measures <- tibble(SSE = c(sum((f_log - f_dc) ^ 2),
               sum((f_log - f_ml) ^ 2)),
       Proc = c(vegan::procrustes((f_log), f_dc)$ss,
                vegan::procrustes((f_log), f_ml)$ss),
       rbind(dr_quality(f_log, f_dc, K = 20)$quality,
             dr_quality(f_log, f_ml, K = 20)$quality)
       )
kable(measures, booktabs = TRUE, digits = 5, escape = FALSE) %>%
  kable_styling(latex_options = "scale_down") %>%
  kable_paper(full_width = FALSE) 
# f_dc gives a better embedding



## --------------------------------------------------------
# Matrix multiplication to compute the Riemannian matrix only around p=c(0,0,1)
## --------------------------------------------------------
# g_q = h_p ^ {-.5} %*% h_q %*% h_p ^ {-.5}
hp <- expm::sqrtm(solve(rmetric[,,N])) # h_N^(-.5)
# all.equal(hp %*% hp, solve(rmetric[,,N]), tolerance = 1e-10)

# r <- 1
fn_dc <- matrix(NA, N, 2)
theta <- rep(NA, N)
for (i in cap_index) {
  gq <- hp %*% rmetric[,,i] %*% hp
  # bindex <- which(adj_matrix[i,] <= r)
  # theta[i] <- sqrt(rmetric_det[i])
  # theta[i] <- sqrt(rmetric_det[i]) ^ (-1)
  # theta[i] <- sqrt(mean(rmetric_det[bindex]) * rmetric_det[i] ) ^ (-1)
  # theta[i] <- sqrt(mean(rmetric_det[bindex]) / rmetric_det[i] ) ^ (-1)
  theta[i] <- sqrt( det(gq) ) # used in the paper, (\frac{|\det \pmb{H}(\pmb{y}_i)|}{|\det \pmb{H}(\pmb{p})|} )^{1/2}
  fn_dc[i,] <- theta[i] %*% fn[i,]
}

# f_metric is not improving
f_log <- rbind(y, c(0,0))[cap_index,]
f_ml <- fn[cap_index,]
f_dc <- fn_dc[cap_index,]
plot(f_log)
points(f_ml, col = "red")
points(f_dc, col = "green")

tibble(SSE = c(sum((f_log - f_dc) ^ 2),
                           sum((f_log - f_ml) ^ 2)),
                   Proc = c(vegan::procrustes((f_log), f_dc)$ss,
                            vegan::procrustes((f_log), f_ml)$ss),
                   rbind(dr_quality(f_log, f_dc, K = 20)$quality,
                         dr_quality(f_log, f_ml, K = 20)$quality)
) %>% 
kable(booktabs = TRUE, digits = 5, escape = FALSE) %>%
  kable_styling(latex_options = "scale_down") %>%
  kable_paper(full_width = FALSE) 
# f_dc does NOT give a better embedding


## --------------------------------------------------------
# Density estimate for the spherical cap and half sphere area
## --------------------------------------------------------
# fixed bandwidth
fn <- metric_isomap$embedding
E1 <- fn[,1] # rename as Ed to match the aesthetics in plot_ellipse()
E2 <- fn[,2]
prob <- c(1, 50, 95, 99) # c(1, 10, 50, 95, 99) #
p_hdr_isomap <- hdrscatterplot_new(E1, E2, kde.package = "ks", levels = prob, noutliers = noutliers, label = NULL)
p_hdr_isomap_p <- p_hdr_isomap$p +
  plot_ellipse(metric_isomap, add = T, ell.no = 50, ell.size = 0,
               color = blues9[5], fill = blues9[1], alpha = 0.2)
h_hdr_isomap <- p_hdr_isomap$den$den$h # same as ks::Hpi.diag(fn)

Rn <- metric_isomap$rmetric # array
adj_matrix <- metric_isomap$adj_matrix
opt.method <- c("AMISE", "MEAN", "SCALED")[2] # ONLY 2 FOR NOW
riem.scale <- 1 # h_hdr_isomap
r <- .1
r <- sqrt(median(apply(Rn, 3, det)))
fisomap <- vkde(x = fn, h = h_hdr_isomap, vh = Rn, r = r, gridsize = gridsize, eval.points = fn, opt.method = opt.method, riem.scale = riem.scale, adj_matrix = adj_matrix) # 0.923

# # check if vkde with grid estimate is the same as hdrcde::interp.2d
# fixgrid_isomap <- vkde(x = fn, h = NULL, gridsize = gridsize)
# summary(fixgrid_isomap)
# interpden_fix <- hdrcde:::interp.2d(fixgrid_isomap$eval.points[[1]], fixgrid_isomap$eval.points[[2]], fixgrid_isomap$estimate, x0 = E1, y0 = E2)
# all.equal(interpden_fix, fisomap$estimate)

p_isomap <- plot_outlier(x = metric_isomap, gridsize = gridsize, prob = prob, riem.scale = riem.scale, f = fisomap, ell.size = 0)
(p_isomap$p + p_hdr_isomap$p) + coord_fixed() + 
  plot_annotation(title = "Left: variable bandwidth; Right: fixed bandwidth", theme = theme(plot.title = element_text(hjust = 0.5)))

cormethod <- c("pearson", "kendall", "spearman")[3]
cor(preswissroll$den, p_isomap$densities, method = cormethod)
cor(preswissroll$den, p_hdr_isomap$densities, method = cormethod)
mean((preswissroll$den - p_isomap$densities) ^ 2)
mean((preswissroll$den - p_hdr_isomap$densities) ^ 2)




# True density in 3d manifold
den_half_sphere %>% summary()
den_cap %>% summary()

# Estimate true density in 2d with metric embedding of log map
fn <- metric_isomap$embedding

## --------------------------------------------------------
# KDE with fixed bandwidth h, use `fn` to estimate the density of theta
h <- Hpi(fn, binned = TRUE)
denks <- ks::kde(x = fn, h = h, eval.points = fn)
fhat_kde <- denks$estimate
summary(fhat_kde)

## --------------------------------------------------------
# DC-KDE, use `fn` to estimate the density of theta
r <- .09
fhat_dckde <- matrix(0, N, N)
for(i in 1:N){
  Hi <- rmetric[,,i]
  # fi <- rep(0, N)
  bindex <- which( (adj_matrix[i, ]) <= r ) # use only neighbors within radius r
  # bindex <- nn.idx[i, 2:(k+1)] # Or use top k neighbors
  for(j in bindex){
    Hj <- rmetric[,,j]
    fhat_dckde[i, j] <- 1 / sqrt( (Hj) / (Hi) ) * ((2 * pi * r^2) ^ (-1/2)) * 
      exp(-1 / 2 / (r^2) * (fn[i] - fn[j]) ^ 2 / Hj) / (pnorm(1) - pnorm(-1)) # for point i, j=y_i[bindex] are used and saved in row i, rowMeans()
  }
  # fi[bindex] <- 1 / sqrt( (cos(acos(fn[bindex]) - acos(fn[i]))) ) * dnorm(x = acos(fn[bindex]), mean = acos(fn[i]), sd = r) / (pnorm(1) - pnorm(-1)) # kernel function with finite bounds [0,1], with true volume density function and geodesic distance
  # fi[bindex] <- 1 / sqrt( det(rmetric) ) * dnorm(x = acos(fn[bindex]), mean = acos(fn[i]), sd = r) / (pnorm(1) - pnorm(-1)) # use estimated Riemannian metric
  # fhat_dckde[, i] <- fhat_dckde[, i] + fi
}
fhat_dckde <- rowMeans(fhat_dckde)
summary(fhat_dckde)

### Calculate the area under the curve, term C
id <- order(fn) # order(theta)
# sum(diff(theta[id]) * rollmean(dentheta[id], 2)) # for true densities
# sum(diff(theta[id]) * rollmean(fhat_kde[id], 2)) # for kde
(AUC <- sum(diff(fn[id]) * rollmean(fhat_dckde[id], 2))) # for DC-KDE
fhat_dckde <- fhat_dckde / AUC # Correct estimates with AUC
summary(fhat_dckde)
summary(dentheta)
summary(fhat_kde)

# bayestestR::area_under_curve(theta, fhat_dckde, method = "trapezoid")
# AUC <- MESS::auc(theta, fhat_dckde, type = "linear") # use approx()

###########################################################
### 7. Comparison: rank correlation, MSE
###########################################################
# Rank correlation
# Pelletier's estimator (ftheta_kderm) is similar to kde (ftheta_kde)
# larger MSE, smaller mean rank error
rbind(comp.den(dentheta, fhat_kde, paste0("KDE ", "h=", round(h,3))), 
      comp.den(dentheta, fhat_dckde, paste0("DC-KDE ", "r=", round(r,3)))
)
# cor(dentheta, fhat_kde, method = "s")
# cor(dentheta, fhat_dckde, method = "s")
# mean((dentheta - fhat_kde) ^ 2)
# mean((dentheta - fhat_dckde) ^ 2)
# # mean((dentheta - fhat_dckde / AUC) ^ 2)
# mean((rank(dentheta) - rank(fhat_kde)) ^ 2)
# mean((rank(dentheta) - rank(fhat_dckde)) ^ 2)

