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
head(X)
range(X)
label <- as.factor(c(rep(1, 0.99*N), rep(2, 0.01*N)))

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

# distance to (0,0,0,0,0) is the same, = radius = 7
# calculate distance to (0,0,0,0) and color them
label <- as.factor(c(rep(1, 0.99*N), rep(2, 0.01*N)))
dist2center <- sqrt(r^2 - X_new[,5]^2)

# save(X_new, N, dist2center, r, den, file = "paper/figures/fived_animation.rda")
# ## Plot 5-d semi-sphere, run once
library(tourr)
colors <- colourvalues::colour_values(- dist2center, palette = "magma") # generate colors for locations
pchs <- c(rep(16, 0.99*N), rep(17, 0.01*N)) # point shapes for kernels
animate_xy(X_new[,1:5], col = colors, pch = pchs, cex = 0.8, fps=15
           # axes = "bottomleft", 
           )
# "figures/tourr_5d_semisphere.png"

# # library(geozoo)
# # sphere <- sphere.hollow(p = 5)
# # sphere$points <- X_new
# # sphere

# # Save gif as animation
# library(gifski)
# png("frame%03d.png")
# par(ask = FALSE)
# # for(i in 1:10)
# #   plot(rnorm(i * 10), main = i)
# animate_xy(X_new[,1:5], col = colors, pch = pchs, cex = 0.8, fps=15
#            # axes = "bottomleft", 
# )
# dev.off()
# png_files <- sprintf("frame%03d.png", 30:80)
# gif_file <- gifski(png_files)
# unlink(png_files)
# utils::browseURL(gif_file)



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
s <- 4 # embedded in 5-D
k <- 20
method <- "Isomap"
annmethod <- "kdtree"
distance <- "euclidean"
treetype <- "kd"
searchtype <- "radius" # change searchtype for radius search based on `radius`, or KNN search based on `k`
radius <- 10 # the bandwidth parameter, \sqrt(\elsilon), as in algorithm. Note that the radius need to be changed for different datasets, not to increase k

riem.scale <- 20 # tune parameter
gridsize <- 20

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
fixden_isomap <- vkde(x = fn, h = NULL, eval.points = fn, binned = FALSE) # gridsize = gridsize
summary(fixden_isomap$estimate)
fixden_isomap$H

## ----VKDE-------------------------------------------------------------------
Rn <- metric_isomap$rmetric # array
# tictoc::tic()
fisomap <- vkde(x = fn, h = Rn*riem.scale, eval.points = fn)
# tictoc::toc()
summary(fisomap$estimate)

# ----compare with true density----------------------------------------------
par(mfrow=c(1,2))
plot(den, fisomap$estimate, main = paste("Variable bandwidth correlation:", round(cor(den, fisomap$estimate), 3)))
plot(den, fixden_isomap$estimate, main = paste("Fixed bandwidth correlation:", round(cor(den, fixden_isomap$estimate), 3)))
cor(den, fisomap$estimate)
cor(den, fixden_isomap$estimate)






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
fixden_lle <- vkde(x = fn, h = NULL, binned = FALSE, eval.points = fn)
summary(fixden_lle$estimate)
fixden_lle$H

## ----VKDE-------------------------------------------------------------------
# tictoc::tic()
flle <- vkde(x = fn, h = Rn*riem.scale, binned = FALSE, eval.points = fn)
# tictoc::toc()
summary(flle$estimate)

# ----compare with true density----------------------------------------------
par(mfrow=c(1,2))
plot(den, flle$estimate, main = paste("Variable bandwidth correlation:", round(cor(den, flle$estimate), 3)))
plot(den, fixden_lle$estimate, main = paste("Fixed bandwidth correlation:", round(cor(den, fixden_lle$estimate), 3)))
cor(den, flle$estimate)
cor(den, fixden_lle$estimate)






# Laplacian eigenmaps

## ----tsne, message=FALSE, warning=FALSE, eval=TRUE------------------------------
x <- train
method <- "LaplacianEigenmaps"
# perplexity <- 30 # round(k / 3) # 30 by default
# theta <- 0 # for exact tSNE in the C++ tSNE Barnes-Hut implementation
# # tictoc::tic()
metric_le <- metricML(x, s = s, k = k, radius = radius, method = method,
                        # annmethod = annmethod, eps = 0, distance = distance, treetype = treetype,
                        searchtype = searchtype,
                        # perplexity = perplexity, theta = theta, 
                      invert.h = TRUE)
fn <- metric_le$embedding
Rn <- metric_le$rmetric
# E1 <- fn[,1]; E2 <- fn[,2]
# fle <- vkde2d(x = E1, y = E2, h = Rn, binned = FALSE, eval.points = fn)
# fxy_le <- hdrcde:::interp.2d(fle$x, fle$y, fle$z, x0 = E1, y0 = E2)
# # plot_embedding(metric_le)
# # plot_ellipse(metric_le, ell.no = 50)
# # plot_contour(metric_le, n.grid = 20, riem.scale = 1/20)
#
#
# ## ---- echo = F------------------------------------------------------------------
# p_le <- plot_outlier(x = metric_le, n.grid = 20, prob = prob, noutliers = 20, riem.scale = 1/20, f = fle, ell.size = 0)
# p_hdr_le <- hdrscatterplot(E1, E2, kde.package = "ks", noutliers = 20)
# p_hdr_le_p <- p_hdr_le +
#   plot_ellipse(metric_le, ell.no = 50, add = T)
# (p_hdr_le_p + p_le$p) + coord_fixed()
#
# # metric_le$embedding <- preswissroll
# # plot_outlier(x = metric_le, n.grid = 20, prob = prob, noutliers = 20, riem.scale = 1/20, f = fle, ell.size = 0)
## ----Fixed bandwidth density estimates with ks::kde-------------------------
fixden_le <- vkde(x = fn, h = NULL, binned = FALSE, eval.points = fn)
summary(fixden_le$estimate)
fixden_le$H

## ----VKDE-------------------------------------------------------------------
Rn <- metric_le$rmetric # array
tictoc::tic()
fle <- vkde(x = fn, h = Rn*riem.scale, binned = FALSE, eval.points = fn)
tictoc::toc()
summary(fle$estimate)

# ----compare with true density----------------------------------------------
par(mfrow=c(1,2))
plot(den, fle$estimate, main = paste("Variable bandwidth correlation:", round(cor(den, fle$estimate), 3)))
plot(den, fixden_le$estimate, main = paste("Fixed bandwidth correlation:", round(cor(den, fixden_le$estimate), 3)))
cor(den, fle$estimate)
cor(den, fixden_le$estimate)



# # tSNE, embeding dimemsion should be either 1, 2 or 3



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
# fumap <- vkde2d(x = E1, y = E2, h = Rn, binned = FALSE, eval.points = fn)
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
fixden_umap <- vkde(x = fn, h = NULL, binned = FALSE, eval.points = fn)
summary(fixden_umap$estimate)
fixden_umap$H

## ----VKDE-------------------------------------------------------------------
Rn <- metric_umap$rmetric # array
# tictoc::tic()
fumap <- vkde(x = fn, h = Rn*riem.scale, binned = FALSE, eval.points = fn)
# tictoc::toc()
summary(fumap$estimate)

# ----check umap----------------------------------------------
par(mfrow=c(1,2))
plot(den, fumap$estimate, main = paste("Variable bandwidth correlation:", round(cor(den, fumap$estimate), 3)))
plot(den, fixden_umap$estimate, main = paste("Fixed bandwidth correlation:", round(cor(den, fixden_umap$estimate), 3)))
cor(den, fumap$estimate)
cor(den, fixden_umap$estimate)







###--------------------------------------------------------
## compareDensity
###--------------------------------------------------------
# first load all .rda files
fxy <- den
label <- as.factor(c(rep(1, 0.99*N), rep(2, 0.01*N)))
dist2center <- sqrt(r^2 - X_new[,5]^2)
methods <- c("isomap", "lle", "le", "umap")
## scatterplot to compare f_xy for ISOMAP
f <- tibble(fxy = fxy, fxy_vkde = fisomap$estimate, fxy_hdr = fixden_isomap$estimate)
f1 <- tibble(fxy = fxy, fxy_vkde = flle$estimate, fxy_hdr = fixden_lle$estimate)
f2 <- tibble(fxy = fxy, fxy_vkde = fle$estimate, fxy_hdr = fixden_le$estimate)
f3 <- tibble(fxy = fxy, fxy_vkde = fumap$estimate, fxy_hdr = fixden_umap$estimate)

summary(f)
cor(f$fxy_vkde, f$fxy)
cor(f$fxy_hdr, f$fxy)
f=f
pf_vkde <- f %>% 
  ggplot(aes(x = fxy, y = fxy_vkde, col = dist2center, shape = label)) + 
  # scale_color_viridis(option = "A", direction = -1) + 
  geom_point() + 
  labs(x = "", y = "ISOMAP", color = "Distance", shape = "Kernels", title = "Variable bandwidth") +
  scale_y_continuous(limits = c(0, max(f$fxy_vkde)), n.breaks = 6)
pf_hdr <- f %>% 
  ggplot(aes(x = fxy, y = fxy_hdr, col = dist2center, shape = label)) + 
  geom_point() + 
  labs(x = "", y = "", color = "Distance", shape = "Kernels", title = "Fixed bandwidth") +
  scale_y_continuous(limits = c(0, max(f$fxy_hdr)), n.breaks = 6)
p <- pf_vkde + pf_hdr

f=f1
pf_vkde <- f %>% 
  ggplot(aes(x = fxy, y = fxy_vkde, col = dist2center, shape = label)) + 
  # scale_color_viridis(option = "A", direction = -1) +  
  geom_point() + 
  labs(x = "", y = "LLE", color = "Distance", shape = "Kernels", title = "") +
  scale_y_continuous(limits = c(0, max(f$fxy_vkde)), n.breaks = 6)
pf_hdr <- f %>% 
  ggplot(aes(x = fxy, y = fxy_hdr, col = dist2center, shape = label)) + 
  # scale_color_viridis(option = "A", direction = -1) +  
  geom_point() + 
  labs(x = "", y = "", color = "Distance", shape = "Kernels", title = "") +
  scale_y_continuous(limits = c(0, max(f$fxy_hdr)), n.breaks = 6)
p1 <- pf_vkde + pf_hdr

f=f2
pf_vkde <- f %>% 
  ggplot(aes(x = fxy, y = fxy_vkde, col = dist2center, shape = label)) + 
  # scale_color_viridis(option = "A", direction = -1) +  
  geom_point() + 
  labs(x = "True density", y = "Laplacian Eigenmaps", color = "Distance", shape = "Kernels", title = "") +
  scale_y_continuous(limits = c(0, max(f$fxy_vkde)), n.breaks = 6)
pf_hdr <- f %>% 
  ggplot(aes(x = fxy, y = fxy_hdr, col = dist2center, shape = label)) + 
  # scale_color_viridis(option = "A", direction = -1) +  
  geom_point() + 
  labs(x = "True density", y = "", color = "Distance", shape = "Kernels", title = "") +
  scale_y_continuous(limits = c(0, max(f$fxy_hdr)), n.breaks = 6)
p2 <- pf_vkde + pf_hdr

f=f3
pf_vkde <- f %>% 
  ggplot(aes(x = fxy, y = fxy_vkde, col = dist2center, shape = label)) + 
  geom_point() + 
  labs(x = "", y = "UMAP", color = "Distance", shape = "Kernels", title = "") +
  scale_y_continuous(limits = c(0, max(f$fxy_vkde)), n.breaks = 6)
pf_hdr <- f %>% 
  ggplot(aes(x = fxy, y = fxy_hdr, col = dist2center, shape = label)) + 
  geom_point() + 
  labs(x = "", y = "", color = "Distance", shape = "Kernels", title = "") +
  scale_y_continuous(limits = c(0, max(f$fxy_hdr)), n.breaks = 6)
p3 <- pf_vkde + pf_hdr

result <- (p/p1/p2/p3) + plot_layout(guides = 'collect') & 
  theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5, face = "bold")) & 
  scale_color_viridis(option = "A", direction = -1) 
gt <- patchwork::patchworkGrob(result)
gt <- gridExtra::grid.arrange(gt, left = "Estimated density", bottom = "True density")
gt
ggsave(paste0("paper/figures/", "fived", N, "_density_comparison_4ml_riem", format(riem.scale, decimal.mark = "_"), ".png"), gt, width = 8, height = 10, dpi = 300)



# Table for density correlations
fxy <- den
dencor <- function(x) cor(x$estimate, fxy)
# dencor(p_lle)
# dencor(p_hdr_lle)
cors <- cbind(
  c(dencor(fisomap), dencor(fixden_isomap)),
  c(dencor(flle), dencor(fixden_lle)),
  c(dencor(fle), dencor(fixden_le)),
  c(dencor(fumap), dencor(fixden_umap))
) 
rownames(cors) <- c("Variable bandwidth", "Fixed bandwidth")
colnames(cors) <- c("ISOMAP", "LLE", "Laplacian Eigenmaps", "UMAP")
cors %>% 
  kableExtra::kbl(caption = "Correlation between true density and estimated density for four manifold learning embeddings.", booktabs = TRUE, digits = 3, escape = FALSE) %>%
  kable_styling(latex_options = "scale_down") %>%
  kable_paper(full_width = TRUE) 

save(cors, file = paste0("paper/figures/CorrelationTable_", "fived", N, "_4ml_riem", format(riem.scale, decimal.mark = "_"), ".rda"))
