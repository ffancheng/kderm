# 5-d meta data
# 1-4: 99% from N(0,1), 1% from N(0, 5^2)
# 5: sqrt(r^2-x1^2-...-x4^2)
# 6-100: 0
# QR decomposition of another random matrix of dimension 100*N, Q: 100*100
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
library(copula)
library(plotly)
library(mvtnorm)
Jmisc::sourceAll(here::here("R/sources"))
set.seed(1)
N <- 200
mu1 = mu2 = 0
s1 <- 1#/(15^2)
s2 <- 5#/(15^2)
p1 <- 0.99
p2 <- 0.01
x1 <- c(rnorm(0.99*N, mu1, s1), rnorm(0.01*N, mu2, s2))
x2 <- c(rnorm(0.99*N, mu1, s1), rnorm(0.01*N, mu2, s2))
x3 <- c(rnorm(0.99*N, mu1, s1), rnorm(0.01*N, mu2, s2))
x4 <- c(rnorm(0.99*N, mu1, s1), rnorm(0.01*N, mu2, s2))
x <- cbind(x1, x2, x3, x4)
head(x)
range(x)

# GMM density
den_calc <- function(x) {
  p1 * dmvnorm(x, rep(mu1, 4), diag(s1, 4)) + 
    p2 *  dmvnorm(x, rep(mu2, 4), diag(s2, 4))
}
den <- apply(x, 1, den_calc)
length(den)
range(den)
summary(den)
# den <- p1 * dmvnorm(x[i,], rep(mu1, 4), diag(s1, 4)) + 
#   p2 *  dmvnorm(x[i,], rep(mu2, 4), diag(s2, 4))

# semi-sphere radius
# scales::rescale(x, c(0,1))
r <- ceiling(max(sqrt(x1^2 + x2^2 + x3^2 + x4^2)))
r
range(x)

# x_new <- cbind(x/r, 
#            x5 = sqrt(1 - (x1^2 + x2^2 + x3^2 + x4^2)/r^2),
#            matrix(0, N, 95)
#            )
x_new <- cbind(x, 
               x5 = sqrt(r^2 - (x1^2 + x2^2 + x3^2 + x4^2)),
               matrix(0, N, 95)
)
# x0 <- x
head(x_new)
x_new %>% as_tibble() %>% summarise(a = x1^2 + x2^2 + x3^2 + x4^2 + x5^2)

# den1 <- apply(x_new[,1:4], 1, den_calc)
# summary(den1)
# all.equal(den, den1) # TRUE

## QR decomposition
# We split any real matrix A into a product A=QR where Q is a matrix with unit norm orthogonal vectors and R is an upper triangular matrix.
# A is a p*N matrix
# Q a p×p matrix with Q'Q=I
# R a p×N upper triangular matrix
# Q'A = R 

# randomly generate another matrix y (different structure from x) of dimension 100*N
# QR decomposition of y to get the rotation matrix
y <- matrix(runif(N*100, 0, 1), 100, N)
head(y)
dim(y)
QR <- qr(y)
QR$rank
Q <- qr.Q(QR); Q
R <- qr.R(QR); R
# We can reconstruct the matrix y from its decomposition as follows:
all.equal(qr.X(QR, complete = TRUE), y[,1:100])
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

train <- t(Q %*% t(x_new))
dim(train)
any(train==0)
all.equal(train %*% t(train), x_new %*% t(x_new))

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
radius <- 20 # the bandwidth parameter, \sqrt(\elsilon), as in algorithm. Note that the radius need to be changed for different datasets, not to increase k


## ----isomap-------------------------------------------------------------
metric_isomap <- metricML(x, s = s, k = k, radius = radius, method = method, invert.h = TRUE, eps = 0,
                          # annmethod = annmethod, distance = distance, treetype = treetype, 
                          searchtype = searchtype
)
summary(metric_isomap)


## ----ggellipse, include=FALSE, eval=FALSE---------------------------------------
## plot_embedding(metric_isomap) +
##   labs(x = "ISO1", y = "ISO2")
## plot_ellipse(metric_isomap, add = F, n.plot = 50, scale = 20,
##              color = blues9[5], fill = blues9[5], alpha = 0.2)


## -------------------------------------------------------------------------------
# fixed bandwidth
fn <- metric_isomap$embedding
dim(fn)


E1 <- fn[,1] # rename as Ed to match the aesthetics in plot_ellipse()
E2 <- fn[,2]
prob <- c(1, 50, 99)
p_hdr_isomap <- hdrscatterplot(E1, E2, levels = prob, noutliers = 20, label = NULL)
p_hdr_isomap_p <- p_hdr_isomap + 
  plot_ellipse(metric_isomap, add = T, n.plot = 50, scale = 100, 
               color = blues9[5], fill = blues9[5], alpha = 0.2)
p_hdr_isomap




## ----outliers-------------------------------------------------------------------
Rn <- metric_isomap$rmetric # array
n.grid <- 10
fisomap <- vkde2d(x = fn[,1], y = fn[,2], h = Rn, n = n.grid)
fxy_isomap <- hdrcde:::interp.2d(fisomap$x, fisomap$y, fisomap$z, x0 = E1, y0 = E2)
# plot_contour(metric_isomap, n.grid = n.grid, scale = 1/20)


## ----hdroutliers----------------------------------------------------------------
# source(here::here("R/sources/hdrplotting.R"))
p_isomap <- plot_outlier(x = metric_isomap, n.grid = 20, prob = prob, scale = 1/8, f = fisomap, ell_size = 0)


## ----compoutlier, eval = FALSE--------------------------------------------------
(p_hdr_isomap_p + p_isomap$p ) + coord_fixed()


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
E1 <- fn[,1]; E2 <- fn[,2]
flle <- vkde2d(x = E1, y = E2, h = Rn, n = n.grid)
fxy_lle <- hdrcde:::interp.2d(flle$x, flle$y, flle$z, x0 = E1, y0 = E2)
# plot_embedding(metric_lle)
# plot_ellipse(metric_lle, n.plot = 50)
# plot_contour(metric_lle, n.grid = 20, scale = 1/20)


## ---- echo = F------------------------------------------------------------------
p_lle <- plot_outlier(x = metric_lle, n.grid = 20, prob = prob, noutliers = 20, scale = 1/20, f = flle, ell_size = 0)
p_hdr_lle <- hdrscatterplot(E1, E2, kde.package = "ks", noutliers = 20)
p_hdr_lle_p <- p_hdr_lle + 
  plot_ellipse(metric_lle, n.plot = 50, add = T)
(p_hdr_lle_p + p_lle$p) + coord_fixed()



# tSNE

## ----tsne, message=FALSE, warning=FALSE, eval=TRUE------------------------------
x <- train
method <- "tSNE"
perplexity <- 30 # round(k / 3) # 30 by default
theta <- 0 # for exact tSNE in the C++ tSNE Barnes-Hut implementation
# tictoc::tic()
metric_tsne <- metricML(x, s = s, k = k, radius = radius, method = method, 
                        # annmethod = annmethod, eps = 0, distance = distance, treetype = treetype, 
                        searchtype = searchtype, 
                        perplexity = perplexity, theta = theta, invert.h = TRUE)
fn <- metric_tsne$embedding
Rn <- metric_tsne$rmetric
E1 <- fn[,1]; E2 <- fn[,2]
ftsne <- vkde2d(x = E1, y = E2, h = Rn, n = n.grid)
fxy_tsne <- hdrcde:::interp.2d(ftsne$x, ftsne$y, ftsne$z, x0 = E1, y0 = E2)
# plot_embedding(metric_tsne)
# plot_ellipse(metric_tsne, n.plot = 50)
# plot_contour(metric_tsne, n.grid = 20, scale = 1/20)


## ---- echo = F------------------------------------------------------------------
p_tsne <- plot_outlier(x = metric_tsne, n.grid = 20, prob = prob, noutliers = 20, scale = 1/20, f = ftsne, ell_size = 0)
p_hdr_tsne <- hdrscatterplot(E1, E2, kde.package = "ks", noutliers = 20)
p_hdr_tsne_p <- p_hdr_tsne + 
  plot_ellipse(metric_tsne, n.plot = 50, add = T)
(p_hdr_tsne_p + p_tsne$p) + coord_fixed()

# metric_tsne$embedding <- preswissroll
# plot_outlier(x = metric_tsne, n.grid = 20, prob = prob, noutliers = 20, scale = 1/20, f = ftsne, ell_size = 0)

# UMAP

## ----umap, message=FALSE, warning=FALSE, eval=TRUE------------------------------
x <- train
method <- "UMAP"
metric_umap <- metricML(x, s = s, k = k, radius = radius, method = method, 
                        # annmethod = annmethod, eps = 0, distance = distance, treetype = treetype, 
                        searchtype = searchtype, 
                        invert.h = TRUE)



## ---- message=FALSE, eval=TRUE--------------------------------------------------
fumap <- metric_umap$embedding
E1 <- fumap[,1]; E2 <- fumap[,2]
fumap <- vkde2d(x = E1, y = E2, h = Rn, n = n.grid)
fxy_umap <- hdrcde:::interp.2d(fumap$x, fumap$y, fumap$z, x0 = E1, y0 = E2)
# plot_embedding(metric_umap)
# plot_ellipse(metric_umap, n.plot = 50)
# plot_contour(metric_umap, n.grid = 20, scale = 1/20)


## ---- echo = F------------------------------------------------------------------
p_umap <- plot_outlier(x = metric_umap, n.grid = 20, prob = prob, noutliers = 20, scale = 1/20, ell_size = 0)
p_hdr_umap <- hdrscatterplot(E1, E2, kde.package = "ks", noutliers = 20)
p_hdr_umap_p <- p_hdr_umap +
  plot_ellipse(metric_umap, n.plot = 50, add = T)
(p_hdr_umap_p + p_umap$p) + coord_fixed()


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
