# 5-d meta data
# 1-4: 99% from N(0,1), 1% from N(0, 5^2)
# 5: sqrt(r^2-x1^2-...-x4^2)
# 6-100: 0
# QR decomposition of another random matrix of dimension 100*100, Q: 100*100
# Rotate t(x) using Q to get rid of zeros
# N*100 dimension data for manifold learning, embedding dimension s=5
## ----libraries, message=FALSE, echo=TRUE, results='hide'------------------------
library(tidyverse)
library(dimRed)
library(reticulate)
library(here)
library(viridis)
library(hdrcde)
library(igraph)
library(matrixcalc)
library(akima)
library(car)
library(ggforce)
library(ks)
library(patchwork)
library(plotly)
library(mvtnorm)
Jmisc::sourceAll(here::here("R/sources"))
set.seed(1234)
N <- 2000
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
head(X)
range(X)

# GMM density as true meta data density
gmmdensity <- function(x) {
  p1 * dmvnorm(x, rep(mu1, 4), diag(s1, 4)) + 
    p2 *  dmvnorm(x, rep(mu2, 4), diag(s2, 4))
}
den <- apply(X, 1, gmmdensity)
length(den)
range(den)
summary(den)
# den <- p1 * dmvnorm(x[i,], rep(mu1, 4), diag(s1, 4)) + 
#   p2 *  dmvnorm(x[i,], rep(mu2, 4), diag(s2, 4))

# semi-sphere radius
r <- ceiling(max(sqrt(x1^2 + x2^2 + x3^2 + x4^2)))
r # radius = 15 for sigma=5
range(X)

# x_new <- cbind(x/r, 
#            x5 = sqrt(1 - (x1^2 + x2^2 + x3^2 + x4^2)/r^2),
#            matrix(0, N, p - 5)
#            )
X_new <- cbind(X, 
               x5 = sqrt(r^2 - (x1^2 + x2^2 + x3^2 + x4^2)),
               matrix(0, N, p - 5)
)
head(X_new)
X_new %>% as_tibble() %>% summarise(x1^2 + x2^2 + x3^2 + x4^2 + x5^2)

# den1 <- apply(x_new[,1:4], 1, den_calc)
# summary(den1)
# all.equal(den, den1) # TRUE

## Plot 5-d semi-sphere
library(tourr)
# col <- RColorBrewer::brewer.pal(3, "Dark2") # unique(flea$species)
# animate_xy(flea[, 1:6], col = flea$species)

animate_xy(X_new[,1:5], col = viridis(length(den), option = "magma"))
# catogorize density levels and color them accordingly

library(geozoo)
sphere <- sphere.hollow(p = 5)
sphere$points <- X_new
sphere



## QR decomposition
# We split any real matrix A into a product A=QR where Q is a matrix with unit norm orthogonal vectors and R is an upper triangular matrix.
# A is a p*N matrix
# Q a p×p matrix with Q'Q=I
# R a p×N upper triangular matrix
# Q'A = R 

# randomly generate another matrix y (different structure from x) of dimension 100*N
# QR decomposition of y to get the rotation matrix
y <- matrix(runif(p*p, 0, 1), p, p)
head(y)
dim(y)
QR <- qr(y)
QR$rank
Q <- qr.Q(QR); Q
R <- qr.R(QR); R
# We can reconstruct the matrix y from its decomposition as follows:
all.equal(qr.X(QR, complete = TRUE), y)
all.equal(t(Q) %*% y, R)

dim(Q)
dim(R)
all.equal(t(Q) %*% Q, diag(QR$rank), tolerance = 1e-5) # =I
det(Q)
Q[1:10,1:10]

# Use Q as the rotation matrix, p*p, Q'Q=I, Q=Q^(-1), det(Q)=+-1
# Rotate t(x) with Q to get a p*N matrix, then transpose to be the input for manifold learning

# A geometric rotation transforms lines to lines, and preserves ratios of distances between points. From these properties it can be shown that a rotation is a linear transformation of the vectors, and thus can be written in matrix form, Qx. 
# A rotation preserves, not just ratios, but distances themselves.
# x'x = (Qx)'(Qx) = x'Q'Qx = x'Ix = x'x

train <- t(Q %*% t(X_new))
dim(train)
sum(train==0)
all.equal(train %*% t(train), X_new %*% t(X_new))





# ----parameters-----------------------------------------------------------------
# Parameters fixed
x <- train
N <- nrow(x)
s <- 5 # embedded in 5-D
k <- 20
method <- "Isomap"
annmethod <- "kdtree"
distance <- "euclidean"
treetype <- "kd"
searchtype <- "radius" # change searchtype for radius search based on `radius`, or KNN search based on `k`
radius <- 10 # the bandwidth parameter, \sqrt(\elsilon), as in algorithm. Note that the radius need to be changed for different datasets, not to increase k

riem.scale <- 10 # tune parameter
gridsize <- 10

## ----isomap-------------------------------------------------------------
metric_isomap <- metricML(x, s = s, k = k, radius = radius, method = method, invert.h = TRUE, eps = 0,
                          # annmethod = annmethod, distance = distance, treetype = treetype, 
                          searchtype = searchtype
)
summary(metric_isomap)

fn <- metric_isomap$embedding
dim(fn)

# E1 <- fn[,1] # rename as Ed to match the aesthetics in plot_ellipse()
# E2 <- fn[,2]
# prob <- c(1, 10, 50, 95, 99)
# p_hdr_isomap <- hdrscatterplot_new(E1, E2, levels = prob, noutliers = 20, label = NULL)
# p_hdr_isomap
# # p_hdr_isomap_p <- p_hdr_isomap$p + 
# #   plot_ellipse(metric_isomap, add = T, ell.no = 50, ell.size = 100, 
# #                color = blues9[5], fill = blues9[5], alpha = 0.2)
# # p_hdr_isomap
# p_hdr_isomap$densities %>% summary()


## ----Fixed bandwidth density estimates with ks::kde-------------------------
fixden_isomap <- vkde(x = fn, h = NULL, gridsize = gridsize, eval.points = fn)
summary(fixden_isomap$estimate)
fixden_isomap$H

## ----VKDE-------------------------------------------------------------------
Rn <- metric_isomap$rmetric # array
# tictoc::tic()
fisomap <- vkde(x = fn, h = Rn*riem.scale, gridsize = gridsize, eval.points = fn)
# tictoc::toc()
summary(fisomap$estimate)

# ----compare with true density----------------------------------------------
par(mfrow=c(1,2))
plot(den, fisomap$estimate, main = paste("Variable bandwidth correlation:", round(cor(den, fisomap$estimate), 3)))
plot(den, fixden_isomap$estimate, main = paste("Fixed bandwidth correlation:", round(cor(den, fixden_isomap$estimate), 3)))
cor(den, fisomap$estimate)
# # [1] 0.903419
cor(den, fixden_isomap$estimate)
# # [1] 0.844388









# LLE

## ----le, message=FALSE, warning=FALSE, eval=TRUE------------------------------
x <- train
method <- "LLE"
metric_lle <- metricML(x, s = s, k = k, radius = radius, method = method, invert.h = TRUE, eps = 0,
                       # annmethod = annmethod, distance = distance, treetype = treetype, 
                       searchtype = searchtype
)
fn <- metric_lle$embedding
Rn <- metric_lle$rmetric

# # E1 <- fn[,1]; E2 <- fn[,2]
# flle <- vkde(x = fn, h = Rn*riem.scale, gridsize = gridsize, eval.points = fn)
# # fxy_lle <- hdrcde:::interp.2d(flle$x, flle$y, flle$z, x0 = E1, y0 = E2)
# # plot_embedding(metric_lle)
# # plot_ellipse(metric_lle, ell.no = 50)
# # plot_contour(metric_lle, n.grid = 20, riem.scale = 1/20)
# 
# 
# ## ---- echo = F------------------------------------------------------------------
# p_lle <- plot_outlier(x = metric_lle, n.grid = 20, prob = prob, noutliers = 20, riem.scale = 1/20, f = flle, ell.size = 0)
# p_hdr_lle <- hdrscatterplot(E1, E2, kde.package = "ks", noutliers = 20)
# p_hdr_lle_p <- p_hdr_lle + 
#   plot_ellipse(metric_lle, ell.no = 50, add = T)
# (p_hdr_lle_p + p_lle$p) + coord_fixed()


## ----Fixed bandwidth density estimates with ks::kde-------------------------
fixden_lle <- vkde(x = fn, h = NULL, gridsize = gridsize, eval.points = fn)
summary(fixden_lle$estimate)
fixden_lle$H

## ----VKDE-------------------------------------------------------------------
tictoc::tic()
flle <- vkde(x = fn, h = Rn*riem.scale, gridsize = gridsize, eval.points = fn)
tictoc::toc()
summary(flle$estimate)

# ----compare with true density----------------------------------------------
par(mfrow=c(1,2))
plot(den, flle$estimate, main = paste("Variable bandwidth correlation:", round(cor(den, flle$estimate), 3)))
plot(den, fixden_lle$estimate, main = paste("Fixed bandwidth correlation:", round(cor(den, fixden_lle$estimate), 3)))
cor(den, flle$estimate)
# # [1] 0.903419
cor(den, fixden_lle$estimate)
# # [1] 0.844388






# # tSNE, embeding dimemsion should be either 1, 2 or 3
# 
# ## ----tsne, message=FALSE, warning=FALSE, eval=TRUE------------------------------
# x <- train
# method <- "tSNE"
# perplexity <- 30 # round(k / 3) # 30 by default
# theta <- 0 # for exact tSNE in the C++ tSNE Barnes-Hut implementation
# # tictoc::tic()
# metric_tsne <- metricML(x, s = 3, k = k, radius = radius, method = method, 
#                         # annmethod = annmethod, eps = 0, distance = distance, treetype = treetype, 
#                         searchtype = searchtype, 
#                         perplexity = perplexity, theta = theta, invert.h = TRUE)
# fn <- metric_tsne$embedding
# Rn <- metric_tsne$rmetric
# # E1 <- fn[,1]; E2 <- fn[,2]
# # ftsne <- vkde2d(x = E1, y = E2, h = Rn, gridsize = gridsize, eval.points = fn)
# # fxy_tsne <- hdrcde:::interp.2d(ftsne$x, ftsne$y, ftsne$z, x0 = E1, y0 = E2)
# # # plot_embedding(metric_tsne)
# # # plot_ellipse(metric_tsne, ell.no = 50)
# # # plot_contour(metric_tsne, n.grid = 20, riem.scale = 1/20)
# # 
# # 
# # ## ---- echo = F------------------------------------------------------------------
# # p_tsne <- plot_outlier(x = metric_tsne, n.grid = 20, prob = prob, noutliers = 20, riem.scale = 1/20, f = ftsne, ell.size = 0)
# # p_hdr_tsne <- hdrscatterplot(E1, E2, kde.package = "ks", noutliers = 20)
# # p_hdr_tsne_p <- p_hdr_tsne + 
# #   plot_ellipse(metric_tsne, ell.no = 50, add = T)
# # (p_hdr_tsne_p + p_tsne$p) + coord_fixed()
# # 
# # # metric_tsne$embedding <- preswissroll
# # # plot_outlier(x = metric_tsne, n.grid = 20, prob = prob, noutliers = 20, riem.scale = 1/20, f = ftsne, ell.size = 0)
# ## ----Fixed bandwidth density estimates with ks::kde-------------------------
# fixden_tsne <- vkde(x = fn, h = NULL, gridsize = gridsize, eval.points = fn)
# summary(fixden_tsne$estimate)
# fixden_tsne$H
# 
# ## ----VKDE-------------------------------------------------------------------
# Rn <- metric_tsne$rmetric # array
# tictoc::tic()
# ftsne <- vkde(x = fn, h = Rn*riem.scale, gridsize = gridsize, eval.points = fn)
# tictoc::toc()
# summary(ftsne$estimate)
# 
# # ----compare with true density----------------------------------------------
# par(mfrow=c(1,2))
# plot(den, ftsne$estimate, main = paste("Variable bandwidth correlation:", round(cor(den, ftsne$estimate), 3)))
# plot(den, fixden_tsne$estimate, main = paste("Fixed bandwidth correlation:", round(cor(den, fixden_tsne$estimate), 3)))
# cor(den, ftsne$estimate)
# # # [1] 0.903419
# cor(den, fixden_tsne$estimate)
# # # [1] 0.844388



# UMAP

## ----umap, message=FALSE, warning=FALSE, eval=TRUE------------------------------
x <- train
method <- "UMAP"
metric_umap <- metricML(x, s = s, k = k, radius = radius, method = method, 
                        # annmethod = annmethod, eps = 0, distance = distance, treetype = treetype, 
                        searchtype = searchtype, 
                        invert.h = TRUE)



## ---- message=FALSE, eval=TRUE--------------------------------------------------
fn <- metric_umap$embedding
# E1 <- fn[,1]; E2 <- fn[,2]
# fumap <- vkde2d(x = E1, y = E2, h = Rn, gridsize = gridsize, eval.points = fn)
# fxy_umap <- hdrcde:::interp.2d(fumap$x, fumap$y, fumap$z, x0 = E1, y0 = E2)
# plot_embedding(metric_umap)
# plot_ellipse(metric_umap, ell.no = 50)
# plot_contour(metric_umap, n.grid = 20, riem.scale = 1/20)


# ## ---- echo = F------------------------------------------------------------------
# p_umap <- plot_outlier(x = metric_umap, n.grid = 20, prob = prob, noutliers = 20, riem.scale = 1/20, ell.size = 0)
# p_hdr_umap <- hdrscatterplot(E1, E2, kde.package = "ks", noutliers = 20)
# p_hdr_umap_p <- p_hdr_umap +
#   plot_ellipse(metric_umap, ell.no = 50, add = T)
# (p_hdr_umap_p + p_umap$p) + coord_fixed()


## ----Fixed bandwidth density estimates with ks::kde-------------------------
fixden_umap <- vkde(x = fn, h = NULL, gridsize = gridsize, eval.points = fn)
summary(fixden_umap$estimate)
fixden_umap$H

## ----VKDE-------------------------------------------------------------------
Rn <- metric_umap$rmetric # array
tictoc::tic()
fumap <- vkde(x = fn, h = Rn*riem.scale, gridsize = gridsize, eval.points = fn)
tictoc::toc()
summary(fumap$estimate)

# ----compare with true density----------------------------------------------
par(mfrow=c(1,2))
plot(den, fumap$estimate, main = paste("Variable bandwidth correlation:", round(cor(den, fumap$estimate), 3)))
plot(den, fixden_umap$estimate, main = paste("Fixed bandwidth correlation:", round(cor(den, fixden_umap$estimate), 3)))
cor(den, fumap$estimate)
# # [1] 0.903419
cor(den, fixden_umap$estimate)
# # [1] 0.844388



## -------------------------------------------------------------------------------
p_den_isomap <- plot_embedding(metric_isomap) + 
  geom_point(aes(col = preswissroll$den))
p_den_lle <- plot_embedding(metric_lle) + 
  geom_point(aes(col = preswissroll$den)) 
p_den_tsne <- plot_embedding(metric_tsne) + 
  geom_point(aes(col = preswissroll$den)) 
p_den_umap <- plot_embedding(metric_umap) + 
  geom_point(aes(col = preswissroll$den)) 


## FINAL plot to use
(
  (((p_den_isomap + ggtitle("ISOMAP")) | (p_den_lle + ggtitle("LLE")) | (p_den_tsne + ggtitle("t-SNE")) | (p_den_umap + ggtitle("UMAP"))) & scale_color_viridis(option = "inferno") & labs(color = "Density") & coord_fixed()) /
    (p_isomap$p | p_lle$p | p_tsne$p | p_umap$p) /
    (p_hdr_isomap | p_hdr_lle | p_hdr_tsne | p_hdr_umap) +  ## add _p for ellipses
    coord_fixed() +
    #     plot_annotation(subtitle = "Top: True densities from Gaussian mixture model;
    # Middle: Outliers using variable kernel density estimate;
    # Bottom: Outliers from `hdrcde` package with fixed bandwidth;") +
    plot_layout(guides = 'collect') ) &
  labs(x = "", y = "") &
  theme(legend.direction = "vertical")

ggsave(paste0("paper/figures/", mapping, "_outliers_comparison_4ml_3cases.png"), width = 8, height = 6, dpi = 500)









# n <- round(N/4)
# R <- matrix(c(1, 0,
#               0, 1), 
#             nrow = 2, ncol = 2)
# # mu1 <- c(7.5, 7.5)
# # mu2 <- c(7.5, 12.5)
# # mu3 <- c(12.5, 7.5)
# # mu4 <- c(12.5, 12.5)
# mu <- matrix(#~mu1, ~mu2,
#   c(7.5, 7.5,
#     7.5, 12.5,
#     12.5, 7.5,
#     12.5, 12.5), 
#   4, 2, byrow = TRUE)
# co <- NULL
# for(i in 1:4) {
#   # mui <- get(paste0("mu", i))
#   # mui <- as.matrix(mu[i,])
#   a <- mvtnorm::rmvnorm(n, mean = mu[i,], sigma = R)
#   # True density is the mean of densities with all four cores (same R different mus)
#   den <- cbind(mvtnorm::dmvnorm(a, mean = mu[1,], sigma = R),
#                mvtnorm::dmvnorm(a, mean = mu[2,], sigma = R),
#                mvtnorm::dmvnorm(a, mean = mu[3,], sigma = R),
#                mvtnorm::dmvnorm(a, mean = mu[4,], sigma = R)
#   ) %>% rowMeans()
#   a <- cbind(a, den, i)
#   co <- rbind(co, a)
# }
# colnames(co) <- c("x", "y", "den", "label")
# co <- co %>%
#   as_tibble() %>% 
#   mutate(# label = as.factor(label),
#     x = (x - 5)/10, y = (y - 5)/10,
#     den = den * 10 # y = h(x) = (x-5)/10; x = 10*y + 5; dx/dy = 10; f_Y(y) = f_X(x) * 10
#   ) %>% 
#   as.matrix()
