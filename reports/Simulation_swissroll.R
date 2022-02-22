## ----setup, include=FALSE-------------------------------------------------------
# knitr::opts_chunk$set(echo = TRUE)


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
# library(copula)
library(plotly)
Jmisc::sourceAll(here::here("R/sources"))
# set.seed(1)

# data size
N <- 2000
p <- 2

# # ----Swiss roll dataset with a sparse area------------------------
# height <- 30
# x <- (3 * pi / 2) * (1 + 2 * runif(N*1.05, 0, 1))
# y <- height * runif(N*1.05, 0 , 1)
# sr <- cbind(x * cos(x), y, x * sin(x))
# colnames(sr) <- c("x", "y", "z")
# sr <- cbind(as_tibble(sr), col = "grey")
# lowden <- x >= 13 & y >= 18
# sr$col[lowden] <- "red"
# removed <- sample(which(lowden), size = N*0.05, replace = FALSE)
# lowdenpoints <- sr[setdiff(which(lowden), removed), ]
# sr <- sr[-removed, ]
# scatterplot3d::scatterplot3d(sr, color = sr$col)
# # rgl::plot3d(sr, col = sr$col)
# 
# 
# 
# # ----Swiss roll from copula-----------------------------------------
# # ## a 3-dimensional normal copula
# # norm.cop <- normalCopula(0.5, dim = p)
# # u <- rCopula(N, norm.cop)
# # den <- dCopula(u, norm.cop)
# 
# # # 3-d t copula
# # t.cop <- tCopula(c(0.5, 0.3), dim = p, dispstr = "toep",
# #                  df = 2, df.fixed = TRUE)
# # u <- rCopula(N, t.cop)
# # den <- dCopula(u, t.cop)
# # 
# ## a 2-dimensional Clayton copula
# cl3 <- claytonCopula(2, dim = p)
# u <- rCopula(N, cl3)
# # pairs(u)
# claytonden <- dCopula(u, cl3) # true density
# # Swiss roll mapping
# u <- (1.5 * pi) * (1 + 2 * u) # U(1.5pi, 4.5pi)
# a <- u[,1]
# claytonsr <- tibble(x = a * cos(a), y = u[,2], z = a * sin(a))
# 
# # highlight true density of u: red for high densities and black for low ones
# col <- rep("grey", N)
# col[head(rank(claytonden), n=20)] <- "red"
# col[tail(rank(claytonden), n=20)] <- "black"
# plot(u[,1:2], xlab = expression(italic(u)[1]*"'"), ylab = expression(italic(u)[2]*"'"), col = col, cex = 0.4)
# scatterplot3d::scatterplot3d(claytonsr, color = col)
# # plot_ly(data = claytonsr, x = ~ x, y = ~ y, z = ~ z, color = claytonden,
# #         type = "scatter3d", mode = "markers", size = 1, text = paste("density:", claytonden))
# 
# 
# 
# # ----Swiss roll from a Gaussian Mixture Model----------------------------------
# # randomly sampling from a Gaussian Mixture Model with centers/means at (7.5,7.5), (7.5,12.5), (12.5,7.5) and (12.5,12.5). The covariance for each gaussian was the 2x2 identity matrix
# # 400 points in each of the four clusters
# 
# ## read data from online resource http://people.cs.uchicago.edu/~dinoj/manifold/swissroll.html
# # preswissroll <- read_table("http://people.cs.uchicago.edu/~dinoj/manifold/preswissroll.dat", col_names = FALSE)
# # preswissroll_label <- read_table("http://people.cs.uchicago.edu/~dinoj/manifold/preswissroll_labels.dat", col_names = FALSE)$X1 %>% as.factor()
# # swissroll <- read_table("http://people.cs.uchicago.edu/~dinoj/manifold/swissroll.dat", col_names = c("x", "y", "z"))
# # 
# # preswissroll %>% 
# #   add_column(label = preswissroll_label) %>% 
# #   ggplot(aes(X1, X2, col = label)) + 
# #   geom_point()
# 
# ## manually generate multivariate normal random numbers
# n <- round(N/4)
# R <- matrix(c(1, 0,
#               0, 1), 
#             nrow = 2, ncol = 2)
# mu1 <- c(7.5, 7.5)
# mu2 <- c(7.5, 12.5)
# mu3 <- c(12.5, 7.5)
# mu4 <- c(12.5, 12.5)
# # mvtnorm::rmvnorm(n, mean = mu, sigma = R)
# # MASS::mvrnorm(n, mu = mu1, Sigma = R)
# preswissroll <- NULL
# for(i in 1:4) {
#   mui <- get(paste0("mu", i))
#   a <- MASS::mvrnorm(n, mu = mui, Sigma = R)
#   den <- mclust::dmvnorm(a, mean = mui, sigma = R)
#   a <- cbind(a, i, den)
#   preswissroll <- rbind(preswissroll, a)
# }
# colnames(preswissroll) <- c("X1", "X2", "label", "den")
# preswissroll <- preswissroll %>%
#   as_tibble() %>% 
#   mutate(label = as.factor(label)) 
# 
# # Swiss Roll mapping (x,y) -> (x cos x, y, x sin x)
# a <- preswissroll$X1
# swissroll <- tibble(x = a * cos(a), y = preswissroll$X2, z = a * sin(a))



###--------------------------------------------------------
## Start from here! 
###--------------------------------------------------------
set.seed(1234)
mapping <- c("Swiss Roll", "semi-sphere", "Twin Peak", "S Curve")[3]
sr <- mldata(N = 2000, meta = "gaussian", mapping = mapping)
swissroll <- sr$data %>% as.data.frame()
preswissroll <- sr$metadata %>% as_tibble() %>% mutate(den = sr$den, label = c(rep(1:4, each = 500)))
colnames(preswissroll) <- c("X1", "X2", "den", "label")
# plot_ly(data = swissroll, x = ~ x, y = ~ y, z = ~ z, color = preswissroll$label,
#         type = "scatter3d", mode = "markers", size = 1)
plot_ly(data = swissroll, x = ~ x, y = ~ y, z = ~ z, color = sr$den,
        type = "scatter3d", mode = "markers", size = 1, text = paste("density:", preswissroll$den))
scatterplot3d::scatterplot3d(swissroll$x, swissroll$y, swissroll$z, color = preswissroll$label, xlab = "X", ylab = "Y", zlab = "Z")



# plotting (run once for different mappings)
metaplot <- preswissroll %>%
  as_tibble() %>% 
  mutate(label = as.factor(label)) %>% 
  ggplot(aes(X1, X2, col = den)) + 
  geom_point() + 
  scale_color_viridis(option = "B") +
  labs(color = "Density")
metaplot
# ggsave("figures/truedensit_4kernels.png", metaplot, height = 6, width = 8, dpi = 500)



# ----parameters-----------------------------------------------------------------
# train <- readRDS(here::here("data/spdemand_1id336tow_train.rds"))
# train <- readRDS(here::here("data/spdemand_3639id_notow_201length.rds"))

train <- swissroll
# Parameters fixed
x <- train
N <- nrow(x)
s <- 2
k <- 20
method <- "Isomap"
annmethod <- "kdtree"
distance <- "euclidean"
treetype <- "kd"
searchtype <- "radius" # change searchtype for radius search based on `radius`, or KNN search based on `k`
radius <- 5 # the bandwidth parameter, \sqrt(\elsilon), as in algorithm. Note that the radius need to be changed for different datasets, not to increase k

riem.scale <- 0.1
gridsize <- 20

## ---- message=FALSE-------------------------------------------------------------
metric_isomap <- metricML(x, s = s, k = k, radius = radius, method = method, invert.h = TRUE, eps = 0,
                          # annmethod = annmethod, distance = distance, treetype = treetype, 
                          searchtype = searchtype
                          )
# summary(metric_isomap)


# ## ----ggellipse, include=FALSE, eval=FALSE---------------------------------------
# plot_embedding(metric_isomap) +
#   labs(x = "ISO1", y = "ISO2")
# plot_ellipse(metric_isomap, add = F, ell.no = 50, ell.size = .5,
#              color = blues9[5], fill = blues9[5], alpha = 0.2)


## -------------------------------------------------------------------------------
# fixed bandwidth
fn <- metric_isomap$embedding
E1 <- fn[,1] # rename as Ed to match the aesthetics in plot_ellipse()
E2 <- fn[,2]
prob <- c(1, 10, 50, 95, 99) # c(1,50,99)
p_hdr_isomap <- hdrscatterplot_new(E1, E2, levels = prob, noutliers = 20, label = NULL)
p_hdr_isomap_p <- p_hdr_isomap$p +
  plot_ellipse(metric_isomap, add = T, ell.no = 50, ell.size = .5,
             color = blues9[5], fill = blues9[1], alpha = 0.2)
p_hdr_isomap


## ----outliers-------------------------------------------------------------------
Rn <- metric_isomap$rmetric # array
# fisomap <- vkde2d(x = fn[,1], y = fn[,2], h = Rn*riem.scale, gridsize = gridsize) # $x $y $z
# fxy_isomap <- hdrcde:::interp.2d(fisomap$x, fisomap$y, fisomap$z, x0 = E1, y0 = E2)
plot_contour(metric_isomap, gridsize = gridsize, riem.scale = riem.scale)

fisomap <- vkde(x = fn, h = NULL, gridsize = gridsize, eval.points = fn)
# str(fisomap)
# summary(fisomap$estimate)

# check if vkde with grid estimate is the same as hdrcde::interp.2d
# fixgrid_isomap <- vkde(x = fn, h = NULL, gridsize = gridsize)
# summary(fixgrid_isomap)
# interpden_fix <- hdrcde:::interp.2d(fixgrid_isomap$eval.points[[1]], fixgrid_isomap$eval.points[[2]], fixgrid_isomap$estimate, x0 = E1, y0 = E2)
# all.equal(interpden_fix, fisomap$estimate)


## ----hdroutliers----------------------------------------------------------------
# source(here::here("R/sources/hdrplotting.R"))
p_isomap <- plot_outlier(x = metric_isomap, gridsize = gridsize, prob = prob, riem.scale = riem.scale, f = fisomap, ell.size = 0)
# all.equal(fxy_isomap, p_isomap$densities)


## ----compoutlier, eval = FALSE--------------------------------------------------
(p_isomap$p + p_hdr_isomap$p) + coord_fixed() + 
  plot_annotation(title = "Left: variable bandwidth; Right: fixed bandwidth", theme = theme(plot.title = element_text(hjust = 0.5)))

# library(profvis)
# profvis({
#   fisomap <- vkde2d(x = fn[,1], y = fn[,2], h = Rn*riem.scale, gridsize = gridsize) # $x $y $z
#   f1 <- vkde(x = fn, h = Rn*riem.scale, gridsize = gridsize)
# })
# 
# str(fisomap)
# str(f1)
# attributes(f1)
# all.equal(fisomap$z, f1$z, check.attributes = F) # TRUE


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
# flle <- vkde2d(x = E1, y = E2, h = Rn*riem.scale, gridsize = gridsize)
# fxy_lle <- hdrcde:::interp.2d(flle$x, flle$y, flle$z, x0 = E1, y0 = E2)
# plot_embedding(metric_lle)
# plot_ellipse(metric_lle, ell.no = 50)
# plot_contour(metric_lle, gridsize = gridsize, riem.scale = 1/20)

flle <- vkde(x = fn, h = NULL, gridsize = gridsize, eval.points = fn)
## ---- echo = F------------------------------------------------------------------
p_lle <- plot_outlier(x = metric_lle, gridsize = gridsize, prob = prob, noutliers = 20, riem.scale = riem.scale, f = flle, ell.size = 0)
p_hdr_lle <- hdrscatterplot_new(E1, E2, kde.package = "ks", noutliers = 20)
# p_hdr_lle_p <- p_hdr_lle + 
#   plot_ellipse(metric_lle, ell.no = 50, add = T)
(p_lle$p + p_hdr_lle$p) + coord_fixed() + 
  plot_annotation(title = "Left: variable bandwidth; Right: fixed bandwidth", theme = theme(plot.title = element_text(hjust = 0.5)))



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
# ftsne <- vkde2d(x = E1, y = E2, h = Rn*riem.scale, gridsize = gridsize)
# fxy_tsne <- hdrcde:::interp.2d(ftsne$x, ftsne$y, ftsne$z, x0 = E1, y0 = E2)
# plot_embedding(metric_tsne)
# plot_ellipse(metric_tsne, ell.no = 50)
plot_contour(metric_tsne, gridsize = gridsize, riem.scale = 1/20)

ftsne <- vkde(x = fn, h = NULL, gridsize = gridsize, eval.points = fn)

## ---- echo = F------------------------------------------------------------------
p_tsne <- plot_outlier(x = metric_tsne, gridsize = gridsize, prob = prob, noutliers = 20, riem.scale = riem.scale, f = ftsne, ell.size = 0)
p_hdr_tsne <- hdrscatterplot_new(E1, E2, kde.package = "ks", noutliers = 20)
p_hdr_tsne_p <- p_hdr_tsne$p + 
  plot_ellipse(metric_tsne, ell.no = 50, add = T)
(p_tsne$p + p_hdr_tsne$p) + coord_fixed() + 
  plot_annotation(title = "Left: variable bandwidth; Right: fixed bandwidth", theme = theme(plot.title = element_text(hjust = 0.5)))

# metric_tsne$embedding <- preswissroll
# plot_outlier(x = metric_tsne, gridsize = gridsize, prob = prob, noutliers = 20, riem.scale = 1/20, f = ftsne, ell.size = 0)

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
E1 <- fn[,1]; E2 <- fn[,2]
fumap <- vkde2d(x = E1, y = E2, h = Rn*riem.scale, gridsize = gridsize)
fxy_umap <- hdrcde:::interp.2d(fumap$x, fumap$y, fumap$z, x0 = E1, y0 = E2)
# plot_embedding(metric_umap)
# plot_ellipse(metric_umap, ell.no = 50)
# plot_contour(metric_umap, gridsize = gridsize, riem.scale = 1/20)


## ---- echo = F------------------------------------------------------------------
p_umap <- plot_outlier(x = metric_umap, gridsize = gridsize, prob = prob, noutliers = 20, riem.scale = riem.scale, ell.size = 0)
p_hdr_umap <- hdrscatterplot_new(E1, E2, kde.package = "ks", noutliers = 20)
p_hdr_umap_p <- p_hdr_umap +
  plot_ellipse(metric_umap, ell.no = 50, add = T)
(p_umap$p + p_hdr_umap$p) + coord_fixed() + 
  plot_annotation(title = "Left: variable bandwidth; Right: fixed bandwidth", theme = theme(plot.title = element_text(hjust = 0.5)))


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
    (p_hdr_isomap$p | p_hdr_lle$p | p_hdr_tsne$p | p_hdr_umap$p) +  ## add _p for ellipses
  coord_fixed() +
#     plot_annotation(subtitle = "Top: True densities from Gaussian mixture model;
# Middle: Outliers using variable kernel density estimate;
# Bottom: Outliers from `hdrcde` package with fixed bandwidth;") +
  plot_layout(guides = 'collect') ) &
  labs(x = "", y = "") &
  theme(legend.direction = "vertical")

ggsave(paste0("paper/figures/", mapping, "_outliers_comparison_4ml_3cases_riem01.png"), width = 8, height = 6, dpi = 500)



## ----compareoutliers---------------------------------------------------------------------------




p_isomap_all <- plot_outlier(x = metric_isomap, gridsize = gridsize, prob = prob, noutliers = N, riem.scale = 1/20, f = fisomap)
p_tsne_all <- plot_outlier(x = metric_tsne, gridsize = gridsize, prob = prob, noutliers =N, riem.scale = 1/20, f = ftsne)
# p_umap_all <- plot_outlier(x = metric_umap, gridsize = gridsize, prob = prob, noutliers =N, riem.scale = 1/20)


den_rank <- order(preswissroll$den)


head(p_isomap_all$outlier, n = 20)
head(p_tsne_all$outlier, n = 20)
# head(p_umap_all$outlier, n = 20)
outlier_isomap <- order(as.numeric(p_isomap_all$outlier))
outlier_tsne <- order(as.numeric(p_tsne_all$outlier))
# outlier_umap <- order(as.numeric(p_umap_all$outlier))
cor(outlier_isomap, outlier_tsne)
plot(outlier_isomap, den_rank, asp = 1)
abline(0, 1)
# topconfects::rank_rank_plot(outlier_isomap, outlier_tsne, n=50)



## -------------------------------------------------------------------------------
head(p_isomap$outlier, n = 20)
head(p_hdr_tsne$outlier, n = 20)
hdroutlier_isomap <- order(as.numeric(p_isomap$outlier))
hdroutlier_tsne <- order(as.numeric(p_hdr_tsne$outlier))
cor(hdroutlier_isomap, hdroutlier_tsne)
plot(hdroutlier_isomap, hdroutlier_tsne, asp = 1)
abline(0, 1)
# cor(hdroutlier_isomap %>% head(50), hdroutlier_tsne %>% head(50))
# topconfects::rank_rank_plot(hdroutlier_isomap, hdroutlier_tsne, n=50)


## -------------------------------------------------------------------------------
tail(p_isomap_all$outlier, n = 20)
tail(p_tsne_all$outlier, n = 20)
cor(outlier_isomap %>% tail(50), outlier_tsne %>% tail(50))


## -------------------------------------------------------------------------------
tail(p_isomap$outlier, n = 20)
tail(p_hdr_tsne$outlier, n = 20)
cor(hdroutlier_isomap %>% tail(50), hdroutlier_tsne %>% tail(50))

