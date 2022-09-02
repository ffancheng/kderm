# rmarkdown::render("R/simulation_s1.R")

rm(list= ls())
library(tidyverse)
library(dimRed)
library(ks)
library(zoo)
Jmisc::sourceAll(here::here("R/sources"))
set.seed(1)

###########################################################
### Data generation
###########################################################
N <- 1000
d <- 1
distribution <- c("uniform", "normal", "beta", "lognormal", "mixed")[4] # CHANGE DISTRIBUTION OF THETA
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
         # 3. theta ~ Beta(a, b) \in [0,1]; (2,5), (5,1)
         # shape1 <- .5
         # shape2 <- .5
         shape1 <- 2
         shape2 <- 5
         theta <- rbeta(N, shape1, shape2)
         dentheta <- dbeta(theta, shape1, shape2)
         theta <- theta * (pi / 2)
         dentheta <- dentheta / (pi / 2)
         # plot(theta, dentheta)
       },
       "lognormal" = {
         meanlog <- 0
         sdlog <- .15
         theta <- rlnorm(N, meanlog, sdlog)
         dentheta <- dlnorm(theta, meanlog, sdlog)
       },
       "mixed" = {
         components <- sample(1:3, prob = c(0.3, 0.5, 0.2), size = 5*N, replace = TRUE)
         mus <- c(.1, .5, 1)
         sds <- sqrt(c(.3, 1, .1) / pi)
         theta <- rnorm(n=N, mean=mus[components], sd=sds[components])
         dentheta <- .3 * dnorm(theta, mean=mus[1], sd=sds[1]) + .5 * dnorm(theta, mean=mus[2], sd=sds[2]) + .2 * dnorm(theta, mean=mus[3], sd=sds[3])
         keep <- theta >= 0 & theta <= pi / 2 # only keep points within 0 and pi/2
         theta <- theta[keep]
         dentheta <- dentheta[keep]
         N <- length(theta)
       })
## 1-D manifold, (theta) \in [0, \pi/4] with dentheta
plot(theta, dentheta)
summary(dentheta)

## 2-D S1 curve, (x, y) with denx
# Transformation
x <- cos(theta) # |dtheta/dx| = |-1/sqrt(1-x^2)|
denx <- dentheta * 1 / sqrt(1 - x ^ 2) # True density of x, f(x) = f(theta) * |dtheta/dx|
y <- sin(theta)
## AUC, check integral to 1
idx <- order(x)
sum(diff(x[idx]) * rollmean(denx[idx],2))
# bayestestR::area_under_curve(x, y, method = "trapezoid")
# MESS::auc(x, y, type = "linear") # use approx()
summary(denx)

###########################################################
### 1. KDE with fixed bandwidth h, use `theta` to estimate the density of theta
###########################################################
ftheta_kde <- kde(theta, eval.points = theta)$estimate
summary(ftheta_kde)
# # This is the same as ftheta_kderm below
h_kde <- hpi(theta, binned = TRUE)
# ftheta_kderm <- rep(0, N)
# for(i in 1:N){
#   fi <- dnorm(x = theta, mean = theta[i], sd = h_kde)
#   ftheta_kderm <- ftheta_kderm + fi
# }
# ftheta_kderm <- ftheta_kderm / N
# summary(ftheta_kderm)

###########################################################
### 2. KDE on Riemannian manifold, use `theta` to estimate the density of theta
###########################################################
# Test if Pelletier's estimator is close to KDE, r is the tuning parameter
r <- .1
ftheta_kderm <- matrix(0, N, N)
for(i in 1:N){
  fi <- rep(0, N)
  bindex <- which( (abs(theta - theta[i]) / r) <= 1 )
  fi[bindex] <- dnorm(x = theta[i], mean = theta[bindex], sd = r) / (pnorm(1) - pnorm(-1)) # for point i, estimate using its NNs with bindex, non-NNs are 0
  ftheta_kderm[i, ] <- ftheta_kderm[i, ] + fi
}
ftheta_kderm <- rowMeans(ftheta_kderm)
summary(ftheta_kderm)

###########################################################
### 3. DC-KDE, use `x` to estimate the density of theta
###########################################################
r <- .08
ftheta_dckde <- matrix(0, N, N)
for(i in 1:N){
  fi <- rep(0, N)
  bindex <- which( (abs(acos(x) - acos(x[i])) / r) <= 1 ) # use only neighbors within radius r
  fi[bindex] <- 1 / sqrt( abs(cos(acos(x[bindex]) - acos(x[i]))) ) * dnorm(x = acos(x[i]), mean = acos(x[bindex]), sd = r) / (pnorm(1) - pnorm(-1)) # kernel function with finite bounds [0,1]
  ftheta_dckde[, i] <- ftheta_dckde[, i] + fi
}
ftheta_dckde <- rowMeans(ftheta_dckde)
summary(ftheta_dckde)

### Calculate the area under the curve, term C
id <- order(acos(x)) # order(theta)
sum(diff(theta[id]) * rollmean(dentheta[id], 2)) # for true densities
sum(diff(theta[id]) * rollmean(ftheta_kde[id], 2)) # for kde
sum(diff(theta[id]) * rollmean(ftheta_kderm[id], 2)) # for Pelletier's estimator
(AUC <- sum(diff(theta[id]) * rollmean(ftheta_dckde[id], 2))) # for DC-KDE
ftheta_dckde <- ftheta_dckde / AUC # Correct estimates with AUC

###########################################################
### 4. Comparison: rank correlation, MSE
###########################################################
# Rank correlation
# Pelletier's estimator (ftheta_kderm) is similar to kde (ftheta_kde)
# larger MSE, smaller mean rank error
comp.den <- function(x, y, row.name = "", cor.method = "s") {
  if (length(x) != length(y)) stop("Length of x and y not matching!")
  mse <- mean((x - y) ^ 2)
  rankcor <- cor(x, y, method = cor.method)
  rankmse <- mean((rank(x) - rank(y)) ^ 2)
  dat <- data.frame(mse = mse, rankcor = rankcor, rankmse = rankmse)
  rownames(dat) <- row.name
  return(dat)
}
rbind(comp.den(dentheta, ftheta_kde, paste0("KDE ", "h=", round(h_kde,3))), 
      comp.den(dentheta, ftheta_dckde, paste0("DC-KDE ", "r=", round(r,3)))
)
# cor(dentheta, ftheta_kde, method = "s")
# cor(dentheta, ftheta_dckde, method = "s")
# mean((dentheta - ftheta_kde) ^ 2)
# mean((dentheta - ftheta_dckde) ^ 2)
# # mean((dentheta - ftheta_dckde / AUC) ^ 2)
# mean((rank(dentheta) - rank(ftheta_kde)) ^ 2)
# mean((rank(dentheta) - rank(ftheta_dckde)) ^ 2)

## Kolmogorov-Smirnov Tests: if both samples are drawn from the same continuous distribution
ks.test(theta, x)
ks.test(dentheta, denx)
ks.test(dentheta, ftheta_kde)
ks.test(dentheta, ftheta_dckde)

## Earth-moving distance between distributions
?emdist
e# move::emd()

## Q-Q plot


###########################################################
### 5. Plotting
###########################################################
par(mfrow=c(1,1))
plot(theta, ftheta_kderm, main = paste("True/Estimated density of theta"), cex = .2, col = "red", lty = 3, xlim = c(min(theta), max(theta)),
     ylim = c(0, max(dentheta, ftheta_kde, ftheta_kderm))
)
points(theta, dentheta, lty = 1, cex = .2, col = 1)
points(theta, ftheta_kde, lty = 2, cex = .2, col = 3)
legend(x = "topright",             # Position
       legend = c("True density", paste("KDE", "h=", round(h_kde,3)), paste("KDERiem", "r=", round(r,3))),  # Legend texts
       lty = c(1, 2, 3),           # Line types
       col = c(1, 3, 2),           # Line colors
       lwd = 2)                    # Line width

par(mfrow=c(2,2))
plot(rank(dentheta), rank(ftheta_kde), main = paste("Rank correlation:", round(cor(rank(dentheta), rank(ftheta_kde), method = "s"), 3)))
plot(rank(dentheta), rank(ftheta_kderm), main = paste("Rank correlation:", round(cor(rank(dentheta), rank(ftheta_kderm), method = "s"), 3)))
plot(dentheta, ftheta_kde); abline(0, 1, lty = "dashed")
plot(dentheta, ftheta_kderm); abline(0, 1, lty = "dashed")

## --------------------------------------------------------
## Summary
# DC-KDE gives estimates with higher rank correlations, but slightly lower MSE compared to KDE up to 1e-4
# The bandwidth scalar r is to be optimized
# The *injectivity radius* essentially measures the size of the largest ball around p for which radial geodesics behave like geodesics in R^d
# Equivalently, it defines the largest ball on which geodesic normal coordinates around p may be used.
# In manifold learning, it should be less than the radius for NN search, injM <= radius


## Optimize r
rs <- seq(0, 1, .01)[-1]
ftheta <- matrix(NA, length(rs), N + 3)
colnames(ftheta) <- c(paste0("f", 1:N), "mse", "rankcor", "rankmse")
for(m in 1:length(rs)){
  r <- rs[m]
  ftheta_dckde <- matrix(0, N, N)
  for(i in 1:N){
    fi <- rep(0, N)
    bindex <- which( (abs(acos(x) - acos(x[i])) / r) <= 1 ) # use only neighbors within radius r
    fi[bindex] <- 1 / sqrt( abs(cos(acos(x[bindex]) - acos(x[i]))) ) * dnorm(x = acos(x[i]), mean = acos(x[bindex]), sd = r) / (pnorm(1) - pnorm(-1)) # kernel function with finite bounds [0,1]
    ftheta_dckde[, i] <- ftheta_dckde[, i] + fi
  }
  ftheta_dckde <- rowMeans(ftheta_dckde)
  ### Calculate the area under the curve, term C
  id <- order(theta)
  (AUC <- sum(diff(theta[id]) * rollmean(ftheta_dckde[id], 2))) # for DC-KDE
  ftheta_dckde <- ftheta_dckde / AUC # Correct estimates with AUC
  ftheta[m, 1:N] <- ftheta_dckde
  ftheta[m, N+1] <- mean((dentheta - ftheta_dckde) ^ 2)
  ftheta[m, N+2] <- cor(dentheta, ftheta_dckde, method = "s")
  ftheta[m, N+3] <- mean((rank(dentheta) - rank(ftheta_dckde)) ^ 2)
}
par(mfrow = c(1, 2))
plot(ftheta[, N + 1])
plot(ftheta[, N + 2])
# plot(ftheta[, N+3])
rs[which.min(ftheta[, "mse"])]
rs[which.max(ftheta[, "rankcor"])]
(r_opt <- rs[floor( mean(c(which.min(ftheta[, "mse"]), which.max(ftheta[, "rankcor"]))) )])

## --------------------------------------------------------
## Optimization rule
# Select r in between the optimized MSE and rank correlation
# floor( mean(which.min(ftheta[,"mse"]), which.min(ftheta[,"rankcor"])) )

## -----------------------------------------------------------------------------------------
###########################################################
### 6. DC-KDE, use ML embedding `fn` from `x` to estimate the density of theta
###########################################################
# Now reduce the dimension and estimate the density of the embedding
train <- data.frame(x = x, y = y)
N <- nrow(train)
s <- 1
k <- 10
method <- "annIsomap"
annmethod <- "kdtree"
distance <- "euclidean"
treetype <- "kd"
searchtype <- "radius" # change searchtype for radius search based on `radius`, or KNN search based on `k`
radius <- 10 # the bandwidth parameter, \sqrt(\elsilon), as in algorithm. Note that the radius need to be changed for different datasets, not to increase k; riemann_metric() requires eigen() on each H
gridsize <- 20
noutliers <- 20
# riem.scale <- .1 # scale Riemmanian

# ISOMAP
suppressWarnings(metric_isomap <- metricML(x = train, s = s, k = k, radius = radius, method = method, invert.h = TRUE, eps = 0,
                                           # annmethod = annmethod, distance = distance, treetype = treetype, 
                                           searchtype = searchtype
                                           ))
# summary(metric_isomap)
fn <- metric_isomap$embedding
# Estimated Riemannian metric
rmetric <- metric_isomap$rmetric
adj_matrix <- metric_isomap$adj_matrix # pairwise distance matrix
nn.idx <- metric_isomap$nn2res$nn.idx # NN index of size N

# opt.method <- "SCALED"
# riem.scale <- 1
# fisomap <- vkde(x = fn, h = Rn, gridsize = gridsize, eval.points = fn, opt.method = opt.method, riem.scale = riem.scale)


## --------------------------------------------------------
# KDE with fixed bandwidth h, use `fn` to estimate the density of theta
h <- hpi(fn, binned = TRUE)
denks <- ks::kde(x = fn, h = h, eval.points = fn)
fhat_kde <- denks$estimate
# summary(fhat_kde)

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


# # # Kolmogorov-Smirnov Tests: if both samples are drawn from the same continuous distribution
# ks.test(theta, x)
# ks.test(dentheta, denx)
# ks.test(dentheta, fhat_kde)
# ks.test(dentheta, fhat_dckde)
## Q-Q plot 
# qqplot(dentheta, fhat_kde)
# qqplot(dentheta, fhat_dckde)

###########################################################
### 8. Plotting
###########################################################
par(mfrow=c(1,1))
plot(theta, fhat_dckde, main = paste("True/Estimated density of theta"), cex = .2, col = "red", lty = 3, xlim = c(min(theta), max(theta)),
     ylim = c(0, max(dentheta, fhat_kde, fhat_dckde))
)
points(theta, dentheta, lty = 1, cex = .2, col = 1)
points(theta, fhat_kde, lty = 2, cex = .2, col = 3)
legend(x = "topright",             # Position
       legend = c("True density", paste("KDE", "h=", round(h,3)), paste("KDERiem", "r=", round(r,3))),  # Legend texts
       lty = c(1, 2, 3),           # Line types
       col = c(1, 3, 2),           # Line colors
       lwd = 2)                    # Line width

par(mfrow=c(2,2))
plot(rank(dentheta), rank(fhat_kde), main = paste("Rank correlation:", round(cor(rank(dentheta), rank(fhat_kde), method = "s"), 3)))
plot(rank(dentheta), rank(fhat_dckde), main = paste("Rank correlation:", round(cor(rank(dentheta), rank(fhat_dckde), method = "s"), 3)))
plot(dentheta, fhat_kde); abline(0, 1, lty = "dashed")
plot(dentheta, fhat_dckde); abline(0, 1, lty = "dashed")


###########################################################
### 9. Parameter optimization
###########################################################
# Use a series of different hs
# Plot MSE and rank correlation for different h, optimize MSE or rank correlations
rs <- seq(0, radius, .01)[-1]
fhat <- matrix(NA, length(rs), N + 3)
colnames(fhat) <- c(paste0("f", 1:N), "mse", "rankcor", "rankmse")
for(m in 1:length(rs)){
  r <- rs[m]
  fhat_dckde <- matrix(0, N, N)
  for(i in 1:N){
    Hi <- rmetric[,,i]
    bindex <- which( (adj_matrix[i, ]) <= r ) # use only neighbors within radius r
    # bindex <- nn.idx[i, 2:(k+1)] # Or use top k neighbors
    for(j in bindex){
      Hj <- rmetric[,,j]
      fhat_dckde[i, j] <- 1 / sqrt( (Hj) / (Hi) ) * ((2 * pi * r^2) ^ (-1/2)) * 
        exp(-1 / 2 / (r^2) * (fn[i] - fn[j]) ^ 2 / Hj) / (pnorm(1) - pnorm(-1)) # for point i, j=y_i[bindex] are used and saved in row i, rowMeans()
    }
  }
  fhat_dckde <- rowMeans(fhat_dckde)
  ### Calculate the area under the curve, term C
  id <- order(theta)
  (AUC <- sum(diff(theta[id]) * rollmean(fhat_dckde[id], 2))) # for DC-KDE
  fhat_dckde <- fhat_dckde / AUC # Correct estimates with AUC
  fhat[m, 1:N] <- fhat_dckde
  fhat[m, N+1] <- mean((dentheta - fhat_dckde) ^ 2)
  fhat[m, N+2] <- cor(dentheta, fhat_dckde, method = "s")
  fhat[m, N+3] <- mean((rank(dentheta) - rank(fhat_dckde)) ^ 2)
}
par(mfrow = c(1, 2))
plot(fhat[, N + 1])
plot(fhat[, N + 2])
# plot(fhat[, N+3])
rs[which.min(fhat[, "mse"])]
rs[which.max(fhat[, "rankcor"])]
opt_index <- floor( mean(c(which.min(fhat[, "mse"]), which.max(fhat[, "rankcor"]))) )
(r_opt <- rs[opt_index])

## --------------------------------------------------------
## Plotting
par(mfrow=c(1,1))
plot(theta, fhat[opt_index, 1:N], main = paste("Estimated density of fn", "h=", round(r_opt, 3)), cex = .2, col = "red", lty = 3, xlim = c(min(theta), max(theta)),
     ylim = c(0, max(dentheta, fhat_kde, fhat[,1:N], na.rm = T))
)
points(theta, dentheta, lty = 1, cex = .2, col = 1)
# text(x = 0.5, y = 0.5, paste("f(theta) = ", round(dentheta[1], 3)), col = "red")
points(theta, fhat_kde, col = 3, lty = 2, cex = .2)
# points(theta, dentheta, col = "red", lty = "dashed", cex = .2)
legend(x = "topright",          # Position
       legend = c("True density", "Estimates with riemannian", "Estimates with kde"),  # Legend texts
       lty = c(1, 2, 3),           # Line types
       col = c(1, 2, 3),           # Line colors
       lwd = 2)                 # Line width
## --------------------------------------------------------
## Comparison
r_range <- seq(which.min(fhat[, "mse"]), which.max(fhat[, "rankcor"]))
comp_opt <- fhat[r_range, N+1:3, drop = FALSE]
rownames(comp_opt) <- paste0("DC-KDE ", "r=", round(rs[r_range],3))
rbind(comp.den(dentheta, fhat_kde, paste0("KDE ", "h=", round(h,3))), 
      comp.den(dentheta, fhat_dckde, paste0("DC-KDE ", "r=", round(r,3))),
      comp_opt
)
