## ----setup, include=FALSE-------------------------------------------------------
# knitr::opts_chunk$set(echo = TRUE)


## ----libraries, message=FALSE, echo=TRUE, results='hide'------------------------
rm(list= ls())
library(tidyverse)
library(dimRed)
# library(reticulate)
# library(here)
library(viridis)
library(hdrcde)
library(igraph)
library(matrixcalc)
# library(akima)
# library(car)
library(ggforce)
library(ks)
library(patchwork)
# library(copula)
library(plotly)
# library(kableExtra)
Jmisc::sourceAll(here::here("R/sources"))

# data size
N <- 2000
p <- 2

###--------------------------------------------------------
## Start from here! 
###--------------------------------------------------------
set.seed(1234)
mapping <- c("Swiss Roll", "semi-sphere", "Twin Peak", "S Curve")[1] # only run [1] and [3]
sr <- mldata(N = N, meta = "gaussian", mapping = mapping)
# sr <- mldata(N = N, meta = "uniform", mapping = mapping)
swissroll <- sr$data %>% as.data.frame()
preswissroll <- sr$metadata %>% as_tibble() %>% mutate(den = sr$den, label = c(rep(1:4, each = N / 4)))
colnames(preswissroll) <- c("X1", "X2", "den", "label")
plot_ly(data = swissroll, x = ~ x, y = ~ y, z = ~ z, color = sr$den, # colored with 2d density
        type = "scatter3d", mode = "markers", size = 1, text = paste("density:", preswissroll$den))
summary(sr$den) # [0, 1]

# Density of 3d twin peaks using change of coordinates
# P_B(w(u,v)) = P_A(u,v) / sqrt(EG-F^2), where
# EG-F^2 = 1 + (pi * cos(pi*u) * tan(3*v)) ^ 2 + 9 * sin(pi * u) ^ 2 / (cos(3 * v)) ^ 4
# P_A(u,v) = sr$den
u <- preswissroll$X1
v <- preswissroll$X2
den_twinpeaks <- sr$den / sqrt(1 + (pi * cos(pi * u) * tan(3 * v)) ^ 2 + 9 * sin(pi * u) ^ 2 / (cos(3 * v)) ^ 4)
summary(den_twinpeaks)
plot_ly(data = swissroll, x = ~ x, y = ~ y, z = ~ z, color = den_twinpeaks, # colored with 3d density
        type = "scatter3d", mode = "markers", size = 1, text = paste("density:", preswissroll$den))




train <- swissroll
# Parameters fixed
x <- train
y <- sr$metadata # given known embedding with known density
y <- NULL
N <- nrow(x)
s <- 2
k <- 20
method <- "Isomap"
annmethod <- "kdtree"
distance <- "euclidean"
treetype <- "kd"
searchtype <- "radius" # change searchtype for radius search based on `radius`, or KNN search based on `k`
radius <- 8 # the bandwidth parameter, \sqrt(\elsilon), as in algorithm. Note that the radius need to be changed for different datasets, not to increase k

gridsize <- 20
noutliers <- 20

## ---- message=FALSE-------------------------------------------------------------
metric_isomap <- metricML(x, fn = y, s = s, k = k, radius = radius, method = method, invert.h = TRUE, eps = 0,
                          # annmethod = annmethod, distance = distance, treetype = treetype, 
                          searchtype = searchtype
)
# summary(metric_isomap)
# all.equal(y, metric_isomap$embedding, check.attributes = F)
fn <- metric_isomap$embedding
Rn <- metric_isomap$rmetric

# ## ----ggellipse, include=FALSE, eval=FALSE---------------------------------------
# plot_embedding(metric_isomap) +
#   labs(x = "ISO1", y = "ISO2")
# plot_ellipse(metric_isomap, add = F, ell.no = 50, ell.size = .1,
#              color = blues9[5], fill = blues9[5], alpha = 0.2)


# Run once with y <- sr$metadata
## -------------------------------------------------------------------------------
# TRUE density of manifold
## -------------------------------------------------------------------------------
# Transformed from sr$den
den_2dmanifold <- sr$den * (apply(Rn, 3, det)) ^ (.5)
summary(den_2dmanifold)
plotmanifold <- preswissroll %>%
  as_tibble() %>%
  mutate(label = as.factor(label)) %>%
  ggplot(aes(X1, X2, col = den_2dmanifold, shape = label)) +
  geom_point() +
  scale_color_viridis(option = "B") +
  labs(color = "Density", shape = "Kernels")
plotmanifold


# THE FOLLOWING ESTIMATION IS BASED ON THE META DATA
## -------------------------------------------------------------------------------
# fixed bandwidth
E1 <- fn[,1] # rename as Ed to match the aesthetics in plot_ellipse()
E2 <- fn[,2]
prob <- c(1, 50, 95, 99) # c(1, 10, 50, 95, 99) #
p_hdr_isomap <- hdrscatterplot_new(E1, E2, kde.package = "ks", levels = prob, noutliers = noutliers, label = NULL)
p_hdr_isomap_p <- p_hdr_isomap$p +
  plot_ellipse(metric_isomap, add = T, ell.no = 50, ell.size = 0,
               color = blues9[5], fill = blues9[1], alpha = 0.2)
# p_hdr_isomap
h_hdr_isomap <- p_hdr_isomap$den$den$h # same as ks::Hpi.diag(fn)
p_hdr_isomap$densities %>% summary()

## ----outliers-------------------------------------------------------------------
# Rn <- metric_isomap$rmetric # array
# fisomap <- vkde2d(x = fn[,1], y = fn[,2], h = Rn*riem.scale, gridsize = gridsize) # $x $y $z
# fxy_isomap <- hdrcde:::interp.2d(fisomap$x, fisomap$y, fisomap$z, x0 = E1, y0 = E2)
# plot_contour(metric_isomap, gridsize = gridsize, riem.scale = riem.scale) # estimate grid densities with vkde()

adj_matrix <- metric_isomap$adj_matrix
opt.method <- c("AMISE", "MEAN", "SCALED")[3] # 3 ONLY FOR NOW, no scaling for Rn
riem.scale <- 1 # h_hdr_isomap
# r <- .1
r <- sqrt(median(apply(Rn, 3, det)))
fisomap <- vkde(x = fn, h = h_hdr_isomap, vh = Rn, r = r, gridsize = gridsize, eval.points = fn, opt.method = opt.method, riem.scale = riem.scale, adj_matrix = adj_matrix) # 0.923
# fisomap <- vkde(x = fn, h = Rn*riem.scale, gridsize = gridsize, eval.points = fn) # 0.64
# # if scaling the bandwidth
# fisomap_MEAN <- vkde(x = fn, h = Rn, gridsize = gridsize, eval.points = fn, opt.method = "MEAN") # 0.964
# fisomap_AMISE <- vkde(x = fn, h = Rn, gridsize = gridsize, eval.points = fn, opt.method = "AMISE") # 0.941
# h_isomap <- fisomap$H
# str(fisomap)
# summary(fisomap$estimate)
# all.equal(fisomap$estimate, p_isomap$densities)

# # check if vkde with grid estimate is the same as hdrcde::interp.2d
# fixgrid_isomap <- vkde(x = fn, h = NULL, gridsize = gridsize)
# summary(fixgrid_isomap)
# interpden_fix <- hdrcde:::interp.2d(fixgrid_isomap$eval.points[[1]], fixgrid_isomap$eval.points[[2]], fixgrid_isomap$estimate, x0 = E1, y0 = E2)
# all.equal(interpden_fix, fisomap$estimate)

## ----hdroutliers----------------------------------------------------------------
p_isomap <- plot_outlier(x = metric_isomap, gridsize = gridsize, prob = prob, riem.scale = riem.scale, f = fisomap, ell.size = 0)
# all.equal(fxy_isomap, p_isomap$densities)

cormethod <- c("pearson", "kendall", "spearman")[3]
trueden <- den_2dmanifold # sr$den
cor(trueden, p_isomap$densities, method = cormethod)
cor(trueden, p_hdr_isomap$densities, method = cormethod)
mean((trueden - p_isomap$densities) ^ 2)
mean((trueden - p_hdr_isomap$densities) ^ 2)

# Plotting
outliers <- head(order(trueden), noutliers)
p_den_isomap <- plot_embedding(metric_isomap) + 
  geom_point(aes(col = trueden)) + 
  scale_color_viridis(option = "inferno") +
  labs(color = "Density") +
  ggplot2::annotate("text",
                    x = metric_isomap$embedding[outliers, 1] + diff(range(metric_isomap$embedding[,1])) / 50, 
                    y = metric_isomap$embedding[outliers, 2] + diff(range(metric_isomap$embedding[,2])) / 50,
                    label = outliers, col = "blue", cex = 2.5 ) 
(p_den_isomap + p_isomap$p + p_hdr_isomap$p) + 
  coord_fixed() + 
  plot_annotation(title = "Left: True manifold; Middle: DC-KDE; Right: KDE", theme = theme(plot.title = element_text(hjust = 0.5)))
# ggsave("~/Desktop/hdrplot_isomap_metafn_density_compare_thetaform1.png", width = 12, height = 6)
ggsave(paste0("~/Desktop/hdrplot_swissroll_isomap_metafn_density_compare_thetaform2_riem", r, "_", format(if(opt.method == "SCALED") riem.scale, decimal.mark = "_"), ".png"), width = 12, height = 6)

