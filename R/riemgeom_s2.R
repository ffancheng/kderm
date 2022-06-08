#' ---
#' title: "Kernel density estimation for S1 using fixed and variable bandwidths"
#' author: "Fan Cheng"
#' output:
#'   html_document:
#'     keep_tex: true
#' ---
#+ echo=T, warning=F, message=F

rm(list= ls())
library(tidyverse)
library(dimRed)
library(ks)
Jmisc::sourceAll(here::here("R/sources"))

set.seed(1)
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
         theta <- rnorm(n=N,mean=mus[components],sd=sds[components])
         dentheta <- .3 * dnorm(theta, mean=mus[1], sd=sds[1]) + .5 * dnorm(theta, mean=mus[2], sd=sds[2]) + .2 * dnorm(theta, mean=mus[3], sd=sds[3])
         keep <- theta >= 0 & theta <= pi / 2 # only keep points within 0 and pi/2
         theta <- theta[keep]
         dentheta <- dentheta[keep]
         N <- length(theta)
       })
plot(theta, dentheta)
# lines(density(abs(theta)))

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
# plot(x, y, asp = 1)
# range(x) # [0, 1]
# plot(x, denx)
# range(denx)


# KDE with fixed bandwidth
h <- ks::hpi(theta, binned = TRUE) # 0.034
# h <- 0.05
fxkde <- kde(theta, eval.points = theta)$estimate
# plot(theta, fxkde)

# KDE with Riemannian metric
# theta_p(q) = sqrt(|cos(dg(p, q|)
# print(h)
h <- .1
d <- 1
fx <- rep(0, N)
# eval.points <- x
for(i in 1:N){
  fi <- rep(0, N)
  bindex <- which((abs(acos(x) - acos(x[i])) / h) <= 1) # use only neighbors within radius r
  fi[bindex] <- 1 / sqrt((cos(acos(x[bindex]) - acos(x[i])))) * dnorm(x = acos(x[5]), mean = acos(x[i]), sd = h) / (pnorm(1) - pnorm(-1))
  # OR use the true theta and true dg and use all data points for a weighted estimation
  # fi <- 1 / sqrt((cos(acos(x) - acos(x[i])))) * dnorm(x = acos(x), mean = acos(x[i]), sd = h)
  fx <- fx + fi
}
fx <- fx / N
# plot(theta, fx, main = "Estimated f(x) with riemannian metric")

# Plotting
# par(mfrow=c(1,1))
plot(theta, dentheta, main = "True density of theta")
# plot(x, denx, main = "True density of x")
points(theta, fxkde, main = paste("Estimated f(x) with KDE, h = ", round(h, 3)), col = "green")
points(theta, fx, main = "Estimated f(x) with riemannian metric", col = "red")
# points(x, dentheta, col = "red", lty = "dashed")

# plot(x[order(x)][1:N*.6], fx[order(x)][1:N*.6], type = "b", cex = .5, main = "Estimated f(x)[1:.6*N] with riemannian metric")
# abline(h = dentheta, col = "red", lty = "dashed")
# text(x = 0.7, y = 0.55, paste("f(theta) = ", round(dentheta[1], 3)), col = "red")

summary(dentheta)
summary(fxkde)
summary(fx)

# rank correlation
# our estimator (fx) is similar to kde (fxkde)
cor(dentheta, fx, method = "s")
cor(dentheta, fxkde, method = "s")
mean((dentheta - fx) ^ 2)
mean((dentheta - fxkde) ^ 2)
mean((rank(dentheta) - rank(fx)) ^ 2)
mean((rank(dentheta) - rank(fxkde)) ^ 2)


# par(mfrow=c(2,2))
# plot(rank(dentheta), rank(fx), main = paste("Rank correlation:", round(cor(rank(dentheta), rank(fx), method = "s"), 3)))
# plot(rank(dentheta), rank(fxkde), main = paste("Rank correlation:", round(cor(rank(dentheta), rank(fxkde), method = "s"), 3)))
# plot(dentheta, fx); abline(0, 1, lty = "dashed")
# plot(dentheta, fxkde); abline(0, 1, lty = "dashed")


# Summary: Pelletier's estimator gives a similar result as the fixed bandwidth KDE.


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
radius <- .5 # the bandwidth parameter, \sqrt(\elsilon), as in algorithm. Note that the radius need to be changed for different datasets, not to increase k; riemann_metric() requires eigen() on each H

gridsize <- 20
noutliers <- 20
# riem.scale <- .1 # .1 # scale Riemmanian

# ISOMAP
## ---- message=FALSE-------------------------------------------------------------
suppressWarnings(metric_isomap <- metricML(x = train, s = s, k = k, radius = radius, method = method, invert.h = TRUE, eps = 0,
                          # annmethod = annmethod, distance = distance, treetype = treetype, 
                          searchtype = searchtype
))
# summary(metric_isomap)

## -------------------------------------------------------------------------------
# KDE with fixed bandwidth
fn <- metric_isomap$embedding
# fn <- fn+0.45
# E1 <- fn[,1] # rename as Ed to match the aesthetics in plot_ellipse()
# xr <- diff(range(fn, na.rm = TRUE))
# xextend <- 0.15
# xr <- c(min(fn) - xr * xextend, max(fn) + xr * xextend)
h <- hpi(fn, binned = TRUE) # 0.0767
denks <- ks::kde(x = fn, h = h, eval.points = fn)
# str(denks)
ffixed <- denks$estimate
summary(ffixed)
summary(dentheta)
cor(dentheta, ffixed, method = "s")
# plot(fn, ffixed)
# points(theta, dentheta, col = "red", cex = 0.3)
# # plot(x, denx)


# Estimated riem metric
rmetric <- metric_isomap$rmetric
opt.method <- "SCALED"
riem.scale <- 1
# fisomap <- vkde(x = fn, h = Rn, gridsize = gridsize, eval.points = fn, opt.method = opt.method, riem.scale = riem.scale)
adj_matrix <- metric_isomap$adj_matrix # NN distance matrix
dim(adj_matrix)
summary(metric_isomap)
nn.idx <- metric_isomap$nn2res$nn.idx

# ### Bandwidth selection https://bookdown.org/egarpor/NP-EAFIT/dens-bwd.html
# ## Plug-in rule
# # Rule-of-thumb
# h <- bw.nrd0(x = fn)
# h <- bw.nrd(x = fn)
# h <- hpi(fn, binned = TRUE) # same as KernSmooth::dpik(fn)
# ## Cross-validation
# h <- bw.ucv(x = fn)
# h <- bw.bcv(x = fn)
# # Polot estimation from Sheather & Jones (1991) 
# h <- bw.SJ(x = fn)

cor(theta, fn) # 1, meaning that fn is a good embedding of theta. Now we estimate the density of fn as an approximation of density estimates of theta.
# The geodesic distance now is fn[i] - fn[j]
# The estimated Riemannian metric for each point i is a scalar
summary(rmetric) # riem metric are smaller than the fixed bandwidth h, too squiggly; 10 times difference
mean(rmetric)
h

# geodesic distance estimation
# Option 1: true geodesic distance d_g
# Option 2: approximate geodesic distance with inner product adjusted by Riem metric (x-x_i)'h_i^(-1)(x-xi)
# Option 3: approximate geodesic distance with shortest path graph distance d_G
d <- 1
fhat <- rep(0, N)
# h <- h # N=100, h=0.147
print(h)
hidet <- apply(rmetric, 3, det) %>% mean()
# hidetmean <- mean(hidet)
# h <- 1 # N(fn[i], sqrt(hi)) 



# h <- hpi(fn, binned = TRUE)
h <- .2
fhat <- rep(0, N)
for(i in 1:N){
  hi <- rmetric[,,i]
  # fi <- dnorm(x = fn, mean = fn[i], sd = h) # KDE with fixed bandwidth, very close to the above fi
  
  # fi <- abs(cos(acos(fn) - acos(fn[i]))) ^ (-1/2) * dnorm(x = acos(fn), mean = acos(fn[i]), sd = h) # true geodesic distance, true volume density function \theta = \sqrt{abs(cos(d_g))}
  # fi <- abs(cos(acos(fn) - acos(fn[i]))) ^ (-1/2) * 1 / h * (2 * pi) ^ (-1/2) * exp(- 1 / 2 / (h^2)  * (acos(fn) - acos(fn[i]))^2)  # true geodesic distance, true volume density function \theta = \sqrt{abs(cos(d_g))}
  
  # fi <- hi ^ (-1/2) * (1 / h) * (2 * pi) ^ (-1/2) * exp(- 1 / 2 * ((acos(fn) - acos(fn[i])) / h) ^ 2)  # true geodesic distance, estimated volume density function sqrt(det(h*i))
  # fi <- hi ^ (-1/2) * dnorm(x = acos(fn), mean = acos(fn[i]), sd = h) # true geodesic distance, estimated Riem metric
  
  # fi <- hi ^ (-1/2) * (1 / h) * (2 * pi) ^ (-1/2) * exp(- 1 / 2 / (h ^ 2) * (fn - fn[i]) ^ 2 / hi)  # estimated geodesic distance, estimated volume density function sqrt(det(h*i))
  # fi <- dnorm(x = fn, mean = fn[i], sd = h * sqrt(hi)) # approximate geodesic distance with inner product adjusted by Riem metric # set h=.2 to give a less wiggly estimation
  
  fi <- rep(0, N)
  # bindex <- adj_matrix[i,] <= h # only smooth over points within the bandwidth radius
  bindex <- (abs(x - x[i]) / h) <= h
  fi[bindex] <- dnorm(x = fn[bindex], mean = fn[i], sd = h * sqrt(hi))
  
  # fi <- hi ^ (-1/2) * (1 / h) * (2 * pi) ^ (-1/2) * exp(- adj_matrix[,i]^2 / (2 * h ^ 2)) # estimated geodesic distance with shortest path graph distance, estimated Riem metric
  
  # # h <- sweep(hi, 1, hidetmean / hi, "*")
  # hineighbor <- mean(rmetric[ , , nn.idx[i, (2:(k+1))]]) # very close to the one above
  # # print(c(hi, hineighbor))
  # hi <- hi * hineighbor / hi
  # fi <- dnorm(x = fn, mean = fn[i], sd = h * sqrt(hi)) # approximate geodesic distance, average hi using the hi's of the neighborhood of fn[i], i.e. hi/sum(hi_neighbors)

  fhat <- fhat + fi
}
fhat <- fhat / N
par(mfrow=c(1,1))
plot(fn, fhat, main = paste("Estimated density of fn", "h=", round(h,3)), cex = .2, col = "red", lty = 3, xlim = c(min(fn), max(theta)),
     ylim = c(0, max(dentheta, ffixed, fhat))
)
points(theta, dentheta, lty = 1, cex = .2, col = 1)
# text(x = 0.5, y = 0.5, paste("f(theta) = ", round(dentheta[1], 3)), col = "red")
points(fn, ffixed, col = 3, lty = 2, cex = .2)
legend(x = "topright",          # Position
       legend = c("True density", "Estimates with riemannian", "Estimates with kde"),  # Legend texts
       # lty = c(1, 2, 3),           # Line types
       col = c(1, 2, 3),           # Line colors
       lwd = 2)                 # Line width




## Different volume density function, \theta = sqrt(det(hi) / det(hj))
# truncated normal distribution, suppK = [0,1]
# h <- hpi(fn, binned = TRUE)
# fhat <- rep(0, N)
f <- matrix(0, N, N)
for(i in 1:N){
  hi <- rmetric[,,i]
  # bindex <- which(adj_matrix[i,] <= h) # only smooth over points within the bandwidth radius
  bindex <- which((abs(fn - fn[i]) / sqrt(hi) / h) <= 1)
  hi <- hi / mean(rmetric[,,bindex])
  for(j in bindex) {
    hj <- rmetric[,,j]
    f[i,j] <- (hj / hi) ^ (-1/2) * (1 / h) * (2 * pi) ^ (-1/2) * exp(-1/2 / (h^2) * (fn[i] - fn[j]) ^ 2 / hi) / (pnorm(1) - pnorm(-1)) # suppK = [0,1], scale to integral to 1
    # unifk <- ifelse((fn[i] - fn[j]) ^ 2 / hi / h <= 1, 1, 0) # uniform kernel
    # f[i,j] <- (hi) ^ (-1/2) * (1 / h) * unifk

  }
}
# fhat <- fhat / N
fhat <- colMeans(f)
summary(colMeans(f))
summary(rowMeans(f))



par(mfrow=c(1,1))
plot(theta, fhat, main = paste("Estimated density of fn", "h=", round(h,3)), cex = .2, col = "red", lty = 3, xlim = c(min(fn), max(theta)),
     ylim = c(0, max(dentheta, ffixed, fhat))
     )
points(theta, dentheta, lty = 1, cex = .2, col = 1)
# text(x = 0.5, y = 0.5, paste("f(theta) = ", round(dentheta[1], 3)), col = "red")
points(theta, ffixed, col = 3, lty = 2, cex = .2)
# points(theta, dentheta, col = "red", lty = "dashed", cex = .2)
legend(x = "topright",          # Position
       legend = c("True density", "Estimates with riemannian", "Estimates with kde"),  # Legend texts
       # lty = c(1, 2, 3),           # Line types
       col = c(1, 2, 3),           # Line colors
       lwd = 2)                 # Line width
# points(density(fn, kernel = "rectangular"), col = "blue", cex = .2)

# text(x = 0.5, y = 0.5, paste("f(theta) = ", round(dentheta[1], 3)), col = "red")

## extract geodesic distances
dg <- matrix(0, N, N)
truedg <- matrix(0, N,N)
for(i in 1:N) {   
  hi <- rmetric[,,i]
  # bindex <- which(adj_matrix[i,] <= h) # only smooth over points within the bandwidth radius
  #bindex <- which(abs(fn - fn[i]) / sqrt(hi) / h <= 1)
  bindex<-1:N
  # hi <- hi / mean(rmetric[,,bindex])
  for(j in bindex) {
    #hj <- rmetric[,,j]
    dg[i,j] <- sqrt((fn[i] - fn[j]) ^ 2 / (0.5*(hi+rmetric[1,1,j])))
    truedg[i,j] <- abs(theta[i] - theta[j])
  }
}

plot(as.vector(truedg)[as.vector(truedg)<0.00001], as.vector(dg)[as.vector(truedg)<0.00001])
abline(0,1)

mean(as.vector(dg) - as.vector(truedg))
# Our estimates: fhat
# KDE estimates: ffixed
# summary(fhat)
# summary(ffixed)
mean((ffixed - dentheta)^2)
mean((fhat - dentheta)^2)

# # How to compare rank of density estimates
cor(dentheta, ffixed, method = "s")
cor(dentheta, fhat, method = "s")
plot(rank(dentheta), rank(fhat), main = paste("Rank correlation:", round(cor(rank(dentheta), rank(fhat), method = "s"), 3)))
plot(rank(dentheta), rank(ffixed), main = paste("Rank correlation:", round(cor(rank(dentheta), rank(ffixed), method = "s"), 3)))


## Use heavy tailed distribution or mixed normal distribution when KDE fails
## h is optimized for kde fixed
## optimize h using MSE for our estimator



# Use a series of different hs
# Plot MSE and rank correlation for different h, optimize MSE or rank correlations
hs <- seq(0, 1, 0.01)[-1]
hmse <- rep(NA, length(hs))
rcors <- rep(NA, length(hs))
fhat <- matrix(0, length(hs), N)

for (k in 1:length(hs)) {
  h <- hs[k]
  
  f <- matrix(0, N, N)
  for(i in 1:N){
    hi <- rmetric[,,i]
    bindex <- which(adj_matrix[i,] <= h) # only smooth over points within the bandwidth radius
    # fi[bindex] <- dnorm(x = fn[bindex], mean = fn[i], sd = h * sqrt(hi))
    
    for(j in bindex){
      hj <- rmetric[,,j]
      f[i,j] <- (hi / hj) ^ (-1/2) * (1 / h) * (2 * pi) ^ (-1/2) * exp(-1/2 / (h^2) * (fn[i] - fn[j]) ^ 2 / hi) / (pnorm(1) - pnorm(0)) # suppK = [0,1], scale to integral to 1
    }
    
  }
  fhat[k,] <- colMeans(f)
  hmse[k] <- mean((fhat[k,] - dentheta)^2)
  rcors[k] <- cor(dentheta, fhat[k,], method = "s")
}
plot(hs, hmse, type = "b")
abline(h = mean((ffixed - dentheta)^2), col = "red")
plot(hs, rcors, type = "b")
abline(h = cor(dentheta, ffixed, method = "s"), col = "red")
hs[which.min(hmse)] # when MSE is minimized
hs[which.max(rcors)] # when rank correlation is maximized
# save(hs, hmse, rcors, theta, dentheta, ffixed, file = "R/riemgeom_s1_hs.rda")

# optimized MSE
par(mfrow=c(1,1))
plot(fn, fhat[which.min(hmse),], main = paste("Estimated density of fn", "h=", round(hs[which.min(hmse)],3)), cex = .2, col = "red", lty = 3, xlim = c(min(fn), max(theta)),
     ylim = c(0, max(dentheta, ffixed, fhat))
)
points(theta, dentheta, lty = 1, cex = .2, col = 1)
# text(x = 0.5, y = 0.5, paste("f(theta) = ", round(dentheta[1], 3)), col = "red")
points(fn, ffixed, col = 3, lty = 2, cex = .2)
# points(theta, dentheta, col = "red", lty = "dashed", cex = .2)
legend(x = "topright",          # Position
       legend = c("True density", "Estimates with riemannian", "Estimates with kde"),  # Legend texts
       lty = c(1, 2, 3),           # Line types
       col = c(1, 2, 3),           # Line colors
       lwd = 2)                 # Line width

# optimized rank correlation
par(mfrow=c(1,1))
plot(fn, fhat[which.min(rcors),], main = paste("Estimated density of fn", "h=", round(hs[which.min(rcors)],3)), cex = .2, col = "red", lty = 3, xlim = c(min(fn), max(theta)),
     ylim = c(0, max(dentheta, ffixed, fhat))
)
points(theta, dentheta, lty = 1, cex = .2, col = 1)
# text(x = 0.5, y = 0.5, paste("f(theta) = ", round(dentheta[1], 3)), col = "red")
points(fn, ffixed, col = 3, lty = 2, cex = .2)
# points(theta, dentheta, col = "red", lty = "dashed", cex = .2)
legend(x = "topright",          # Position
       legend = c("True density", "Estimates with riemannian", "Estimates with kde"),  # Legend texts
       lty = c(1, 2, 3),           # Line types
       col = c(1, 2, 3),           # Line colors
       lwd = 2)                 # Line width

# Another h value in between
j <- 50
par(mfrow=c(1,1))
plot(fn, fhat[j,], main = paste("Estimated density of fn", "h=", round(hs[j],3)), cex = .2, col = "red", lty = 3, xlim = c(min(fn), max(theta)),
     ylim = c(0, max(dentheta, ffixed, fhat))
)
points(theta, dentheta, lty = 1, cex = .2, col = 1)
# text(x = 0.5, y = 0.5, paste("f(theta) = ", round(dentheta[1], 3)), col = "red")
points(fn, ffixed, col = 3, lty = 2, cex = .2)
# points(theta, dentheta, col = "red", lty = "dashed", cex = .2)
legend(x = "topright",          # Position
       legend = c("True density", "Estimates with riemannian", "Estimates with kde"),  # Legend texts
       lty = c(1, 2, 3),           # Line types
       col = c(1, 2, 3),           # Line colors
       lwd = 2)                 # Line width
