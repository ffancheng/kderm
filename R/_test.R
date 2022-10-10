set.seed(1)
N <- 2000
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
# semi-sphere radius
# scales::rescale(x, c(0,1))
r <- ceiling(max(sqrt(x1^2 + x2^2 + x3^2 + x4^2)))
r
range(X)

# X <- matrix(rnorm(4*N,0,1), N, 4)
X_new <- cbind(X, 
               x5 = sqrt(6^2 - (x1^2 + x2^2 + x3^2 + x4^2))
)
X_new %>% as_tibble() %>% summarise(x1^2 + x2^2 + x3^2 + x4^2 + x5^2) # radius is 15

library(geozoo)
sphere <- sphere.hollow(p = 5)
sphere$points <- X_new
sphere

# GMM density as true meta data density
gmmdensity <- function(x) {
  p1 * dmvnorm(x, rep(mu1, 4), diag(s1, 4)) + 
    p2 *  dmvnorm(x, rep(mu2, 4), diag(s2, 4))
}
den <- apply(X, 1, gmmdensity)
length(den)
range(den)
summary(den)

library(tourr)
animate_xy(X_new[,1:5], col = viridis(length(den),option = "magma"))

col <- RColorBrewer::brewer.pal(length(den), "Dark2") # unique(flea$species)
animate_xy(flea[, 1:6], col = flea$species)




library(ks)
set.seed(8192)
samp <- 200
mus <- rbind(c(-2,2), c(0,0), c(2,-2))
Sigmas <- rbind(diag(2), matrix(c(0.8, -0.72, -0.72, 0.8), nrow=2), diag(2))
cwt <- 3/11
props <- c((1-cwt)/2, cwt, (1-cwt)/2)
x <- rmvnorm.mixt(n=samp, mus=mus, Sigmas=Sigmas, props=props)

plotmixt(mus=mus, Sigmas=Sigmas, props=props, xlim=c(-4,4), ylim=c(-4,4))
plot(x, xlim=c(-4,4), ylim=c(-4,4), xlab="x", ylab="y")








set.seed(1)
N <- 200
n <- round(N/4)
R <- matrix(c(.02, 0,
              0, .02), # variance-covariance matrix
            nrow = 2, ncol = 2)
# mu <- matrix(#~mu1, ~mu2,
#               c(7.5, 7.5,
#               7.5, 12.5,
#               12.5, 7.5,
#               12.5, 12.5), 
#               4, 2, byrow = TRUE)
mu <- matrix(#~mu1, ~mu2,
  c(0.25, 0.25,
    0.25, 0.75,
    0.75, 0.25,
    0.75, 0.75), 
  4, 2, byrow = TRUE)
co <- NULL
for(i in 1:4) {
  a <- mvtnorm::rmvnorm(n, mean = mu[i,], sigma = R)
  # True density is the mean of densities with all four cores (same R different mus)
  den <- cbind(mvtnorm::dmvnorm(a, mean = mu[1,], sigma = R),
               mvtnorm::dmvnorm(a, mean = mu[2,], sigma = R),
               mvtnorm::dmvnorm(a, mean = mu[3,], sigma = R),
               mvtnorm::dmvnorm(a, mean = mu[4,], sigma = R)
  ) %>% rowMeans()
  a <- cbind(a, den, i)
  co <- rbind(co, a)
}
colnames(co) <- c("x", "y", "den", "label")

set.seed(1)
co <- NULL
for(i in 1:4) {
  a <- mvtnorm::rmvnorm(n, mean = mu[i,], sigma = R)
  # # True density is the mean of densities with all four cores (same R different mus)
  den <- cbind(mvtnorm::dmvnorm(a, mean = mu[1,], sigma = R),
               mvtnorm::dmvnorm(a, mean = mu[2,], sigma = R),
               mvtnorm::dmvnorm(a, mean = mu[3,], sigma = R),
               mvtnorm::dmvnorm(a, mean = mu[4,], sigma = R)
  ) %>% rowMeans()
  # den <- ks::dmvnorm.mixt(a, mus = mu, Sigmas = rbind(R,R,R,R), props = rep(1,4)/4)
  a <- cbind(a, den, i)
  co <- rbind(co, a)
}
colnames(co) <- c("x", "y", "den", "label")


# check area=1
E1 <- co[,1]; E2 <- co[,2]
f <- hdrcde:::den.estimate.2d(x = E1, y = E2, kde.package = "ks", xextend=0.15, yextend = 0.15)
str(f)
range(f$z)
mean(diff(f$x) * diff(f$y)) * sum(f$z) # 0.965

# all code needed to check AUC=1
a <- mldata(200, 2, "gaussian", "Twin Peak")
E1 <- a$metadata[,1]; E2 <- a$metadata[,2]
f <- hdrcde:::den.estimate.2d(x = E1, y = E2, kde.package = "ks", xextend=0.15, yextend = 0.15)
str(f)
range(f$z)
mean(diff(f$x) * diff(f$y)) * sum(f$z) # 0.965
filled.contour(f)
summary(a$den)















a <- p_hdr_isomap$den$den
str(a)
# a$x
# a$y
# a$z
b <- list(x=a$x, y=a$y, z=a$z)
filled.contour(b)
mean(diff(a$x) * diff(a$y)) * sum(a$z) # 0.996

fisomap %>% str()
filled.contour(fisomap)
mean(diff(fisomap$eval.points[[1]]) * diff(fisomap$eval.points[[2]])) * sum(fisomap$estimate) # 0.965


out <- p_isomap$hdr2d_info$den 
filled.contour(out)
mean(diff(out$x) * diff(out$y)) * sum(out$z) # 0.987

summary(a$z)
all.equal(a$z, out$z)

all.equal(fxy_isomap, p_hdr_isomap$densities)

# compare with true density
par(mfrow=c(1,2))
plot(preswissroll$den, p_isomap$densities)
plot(preswissroll$den, p_hdr_isomap$densities)

cor(preswissroll$den, p_isomap$densities)
# [1] 0.9119917
cor(preswissroll$den, p_hdr_isomap$densities)
# [1] 0.6264214

cor(preswissroll$den, fxy_isomap)




## bivariate mixtures 
mus <- rbind(c(-1,0), c(1, 2/sqrt(3)), c(1,-2/sqrt(3)))
Sigmas <- 1/25*rbind(invvech(c(9, 63/10, 49/4)), invvech(c(9,0,49/4)), invvech(c(9,0,49/4)))
props <- c(3,3,1)/7
# dfs <- c(7,3,2)
x <- rmvnorm.mixt(10, mus=mus, Sigmas=Sigmas, props=props)
f <- dmvnorm.mixt(x, mus, Sigmas, props)
all.equal(f, x)

a <- list(x=x[,1], y=x[,2], z=f)
filled.contour(a)
mean(diff(a$x) * diff(a$y)) * sum(a$z) # 0.996


