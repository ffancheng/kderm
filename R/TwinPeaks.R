rm(list = ls())
library(tidyverse)
library(dimRed)
library(viridis)
library(hdrcde)
library(igraph)
library(matrixcalc)
library(ggforce)
library(ks)
library(patchwork)
library(plotly)
library(kableExtra)
Jmisc::sourceAll(here::here("R/sources"))

# Data size
N <- 1000
p <- 2

###--------------------------------------------------------
## Start from here! 
###--------------------------------------------------------
set.seed(1234)
mapping <- c("Swiss Roll", "semi-sphere", "Twin Peak", "S Curve")[3] # only run [1] and [3]
sr <- mldata(N = N, meta = "gaussian", mapping = mapping)
# sr <- mldata(N = N, meta = "uniform", mapping = mapping)
swissroll <- sr$data %>% as.data.frame()
preswissroll <- sr$metadata %>% as_tibble() %>% mutate(den = sr$den, label = c(rep(1:4, each = N / 4)))
colnames(preswissroll) <- c("X1", "X2", "den", "label")
plot_ly(data = swissroll, x = ~ x, y = ~ y, z = ~ z, color = sr$den, # colored with 2d meta data density
        type = "scatter3d", mode = "markers", size = 1, text = paste("density:", preswissroll$den))
summary(sr$den) # [0, 1]

# Density of 3d twin peaks using change of coordinates
# P_B(w(u,v)) = P_A(u,v) / sqrt(EG-F^2), where
# EG-F^2 = 1 + (pi * cos(pi*u) * tan(3*v)) ^ 2 + 9 * sin(pi * u) ^ 2 / (cos(3 * v)) ^ 4
# P_A(u,v) = sr$den
u <- preswissroll$X1
v <- preswissroll$X2
if(mapping == "Swiss Roll") {
  den_twinpeaks <- sr$den / sqrt(1 + u ^ 2)
} else if(mapping == "Twin Peak") {
  den_twinpeaks <- sr$den / sqrt(1 + (pi * cos(pi * u) * tan(3 * v)) ^ 2 + 9 * sin(pi * u) ^ 2 / (cos(3 * v)) ^ 4)
} else {
  print("Double check mldata() for the density after change of coordinates!")
}
summary(den_twinpeaks)
plot_ly(data = swissroll, x = ~ x, y = ~ y, z = ~ z, color = den_twinpeaks, # colored with 3d density
        type = "scatter3d", mode = "markers", size = 1, text = paste("density:", preswissroll$den))




train <- swissroll
# Parameters fixed
x <- train
y <- sr$metadata # given known embedding with known density
# y <- NULL
N <- nrow(x)
s <- 2
k <- N / 20
method <- "Isomap"
annmethod <- "kdtree"
distance <- "euclidean"
treetype <- "kd"
searchtype <- "priority" # change searchtype for radius search based on `radius`, or KNN search based on `k`
radius <- 8 # the bandwidth parameter, \sqrt(\elsilon), as in algorithm. Note that the radius need to be changed for different datasets, not to increase k

gridsize <- 20
noutliers <- 20

metric_meta <- metricML(x, fn = y, s = s, k = k, radius = radius, method = method, invert.h = TRUE, eps = 0,
                          # annmethod = annmethod, distance = distance, treetype = treetype, 
                          searchtype = searchtype
)
# summary(metric_meta)
# all.equal(y, metric_meta$embedding, check.attributes = F)
fn <- metric_meta$embedding
Rn <- metric_meta$rmetric

# ## ----ggellipse, include=FALSE, eval=FALSE---------------------------------------
# plot_embedding(metric_meta) +
#   labs(x = "ISO1", y = "ISO2")
# plot_ellipse(metric_meta, add = F, ell.no = 500, ell.size = .1,
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
  ggplot(aes(X1, X2, col = den_2dmanifold # , shape = label
             )) +
  geom_point() +
  scale_color_viridis(option = "B") +
  labs(color = "Density" # , shape = "Kernels"
       )
plotmanifold
# plot_ly(data = swissroll, x = ~ x, y = ~ y, z = ~ z, color = den_2dmanifold, # colored with 3d density
#         type = "scatter3d", mode = "markers", size = 1, text = paste("density:", preswissroll$den))
trueden <- den_2dmanifold # sr$den


# 0. THE FOLLOWING ESTIMATION IS BASED ON THE META DATA
## -------------------------------------------------------------------------------
# KDE with KDE
E1 <- fn[,1] # rename as Ed to match the aesthetics in plot_ellipse()
E2 <- fn[,2]
prob <- c(1, 50, 95, 99) # c(1, 10, 50, 95, 99) #
p_hdr_meta <- hdrscatterplot_new(E1, E2, kde.package = "ks", levels = prob, noutliers = noutliers, label = NULL)
h_hdr_meta <- p_hdr_meta$den$den$h # same as ks::Hpi.diag(fn)
p_hdr_meta$densities %>% summary()

# DC-KDE
adj_matrix <- metric_meta$adj_matrix
opt.method <- c("AMISE", "MEAN", "SCALED")[3] # 3 ONLY FOR NOW, no scaling for Rn
riem.scale <- 1 # h_hdr_meta
# r <- .1
r <- sqrt(mean(apply(Rn, 3, det)))
fmeta <- vkde(x = fn, h = h_hdr_meta, vh = Rn, r = r, gridsize = gridsize, eval.points = fn, opt.method = opt.method, riem.scale = riem.scale, adj_matrix = adj_matrix)
p_meta <- plot_outlier(x = metric_meta, gridsize = gridsize, prob = prob, riem.scale = riem.scale, f = fmeta, ell.size = 0)


# Rank correlation
cormethod <- c("pearson", "kendall", "spearman")[3]
trueden <- den_2dmanifold # sr$den
cor(trueden, p_meta$densities, method = cormethod)
cor(trueden, p_hdr_meta$densities, method = cormethod)
mean((trueden - p_meta$densities) ^ 2)
mean((trueden - p_hdr_meta$densities) ^ 2)
# truerank <- rank(trueden) # top 10 are the 10 anomalies
# rank_dckde <- rank(p_meta$densities)
# rank_kde <- rank(p_hdr_meta$densities)
# # cor(truerank, rank_dckde)

f <- tibble(fxy = trueden, # true densities
            fxy_vkde = p_meta$densities, fxy_hdr = p_hdr_meta$densities
)
# # Plot ranks instead of densities
# uncomment when plotting ranks instead of densities
f <- f %>% summarise_all(rank)
# cor(f$fxy_vkde, f$fxy, method = "spearman")
# cor(f$fxy_hdr, f$fxy, method = "spearman")


pf_vkde <- f %>% 
  ggplot(aes(x = fxy, y = fxy_vkde, col = factor(preswissroll$label), shape = factor(preswissroll$label))) + 
  geom_point() + 
  labs(x = "", y = "", color = "Kernels", shape = "Kernels", title = paste("DC-KDE correlation", round(cor(f$fxy_vkde, f$fxy, method = "spearman"), 3)) ) +
  scale_y_continuous(n.breaks = 6)
pf_hdr <- f %>% 
  ggplot(aes(x = fxy, y = fxy_hdr, col = factor(preswissroll$label), shape = factor(preswissroll$label))) + 
  geom_point() + 
  labs(x = "", y = "", color = "Kernels", shape = "Kernels", title = paste("KDE correlation", round(cor(f$fxy_hdr, f$fxy, method = "spearman"), 3))) +
  scale_y_continuous(limits = c(0, max(f$fxy_hdr)), n.breaks = 5)
# Rank plot
pf_vkde + pf_hdr




# Plotting
trueoutliers <- head(order(trueden), noutliers)
# p_den_meta <- plot_embedding(metric_meta) +
#   geom_point(aes(col = trueden)) +
#   scale_color_viridis(option = "inferno") +
#   labs(color = "Density") +
#   ggplot2::annotate("text",
#                     x = metric_meta$embedding[trueoutliers, 1] + diff(range(metric_meta$embedding[,1])) / 50,
#                     y = metric_meta$embedding[trueoutliers, 2] + diff(range(metric_meta$embedding[,2])) / 50,
#                     label = trueoutliers, col = "blue", cex = 2.5 )

# Categorize density to 5 levels prob
plot_embedding_color_levels <- function(x = trueden, N, prob = c(1, 50, 95, 99), noutliers, fn = NULL) {
  if(max(prob) > 50){ # Assume prob is coverage percentage
    alpha <- (100 - prob) / 100
  } else {# prob is tail probability (for backward compatibility)
    alpha <- prob
  }
  alpha <- sort(alpha)
  falpha <- quantile(x, alpha)
  region <- numeric(N) + 100
  for (i in seq_along(prob)) region[x > falpha[i]] <- 100 - prob[i]
  if (is.null(noutliers)) noutliers <- sum(region > max(prob))
  noutliers <- min(noutliers, N)
  
  # Construct region factor
  prob <- sort(unique(region[region < 100]), decreasing = TRUE)
  prob <- c(prob, 100)
  Region <- factor(region,
                   levels = prob,
                   labels = c(paste(head(prob, -1)), ">99")
  )
  # Sort data so the larger regions go first (other than outliers)
  rr <- region
  rr[region == 100] <- 0
  ord <- order(rr, decreasing = TRUE)
  fn <- fn[ord,] # data is reordered by ord
  
  if (noutliers > 0) {
    outlier_rank <- order(x[ord])  # order den$fxy as well to match the new data
    outliers <- outlier_rank[seq_len(noutliers)] # take top noutliers labels for annotation
  }
  
  p <- fn %>% 
    as_tibble() %>% 
    ggplot(aes(x = E1, y = E2)) + 
    ggplot2::geom_point(ggplot2::aes_string(col = "Region")) + 
    ggplot2::scale_colour_manual(
      name = "HDRs",
      breaks = c(paste(head(sort(prob), -1)), ">99"),
      values = c(RColorBrewer::brewer.pal(length(prob), "YlOrRd")[-1], "#000000")) + 
    ggplot2::annotate("text",
                      x = fn[outliers, 1] + diff(range(fn[outliers, 1])) / 50, 
                      y = fn[outliers, 2] + diff(range(fn[outliers, 2])) / 50,
                      label = seq_len(N)[ord[outliers]], col = "blue", cex = 2.5
    )
  
  return(p)
}

p_den_meta <- plot_embedding_color_levels(x = trueden, N, prob, noutliers, fn = metric_meta$embedding)

(p_den_meta + p_meta$p + p_hdr_meta$p) + 
  coord_fixed() + 
  plot_annotation(title = "Left: True manifold density; Middle: DC-KDE; Right: KDE", theme = theme(plot.title = element_text(hjust = 0.5)))
# ggsave("~/Downloads/hdrplot_meta_metafn_density_compare.png", width = 12, height = 6)
# ggsave(paste0("~/Downloads/hdrplot_swissroll_meta_metafn_density_compare_thetaform2_riem", r, "_", format(if(opt.method == "SCALED") riem.scale, decimal.mark = "_"), ".png"), width = 12, height = 6)

  




# 1. THE FOLLOWING ESTIMATION IS BASED ON THE ISOMAP EMBEDDING
## -------------------------------------------------------------------------------
x <- train
y <- NULL
searchtype <- "radius"
method <- "Isomap"
metric_isomap <- metricML(x, fn = y, s = s, k = k, radius = radius, method = method, invert.h = TRUE, eps = 0,
                          annmethod = annmethod, distance = distance, treetype = treetype,
                          searchtype = searchtype
                          )
fn <- metric_isomap$embedding
Rn <- metric_isomap$rmetric
adj_matrix <- metric_isomap$adj_matrix

## -------------------------------------
# KDE with KDE
E1 <- fn[,1] # rename as Ed to match the aesthetics in plot_ellipse()
E2 <- fn[,2]
# prob <- c(1, 50, 95, 99) # c(1, 10, 50, 95, 99) #
p_hdr_isomap <- hdrscatterplot_new(E1, E2, kde.package = "ks", levels = prob, noutliers = noutliers, label = NULL)
h_hdr_isomap <- p_hdr_isomap$den$den$h # same as ks::Hpi.diag(fn)
p_hdr_isomap$densities %>% summary()

# DC-KDE
# opt.method <- c("AMISE", "MEAN", "SCALED")[3] # 3 ONLY FOR NOW, no scaling for Rn
# riem.scale <- 1 # h_hdr_isomap
# r <- .1
# r <- sqrt(median(apply(Rn, 3, det)))
fisomap <- vkde(x = fn, h = h_hdr_meta, vh = Rn, r = r, gridsize = gridsize, eval.points = fn, opt.method = opt.method, riem.scale = riem.scale, adj_matrix = adj_matrix)
p_isomap <- plot_outlier(x = metric_isomap, gridsize = gridsize, prob = prob, riem.scale = riem.scale, f = fisomap, ell.size = 0)

# plot_ellipse(metric_isomap, add = F, ell.no = 500, ell.size = .1,
#              color = blues9[5], fill = blues9[5], alpha = 0.2)

p_den_isomap <- plot_embedding_color_levels(p = plot_embedding(metric_isomap), x = trueden, N, prob, noutliers, fn = metric_isomap$embedding)
(p_den_isomap + p_isomap$p + p_hdr_isomap$p)+ 
  coord_fixed()


# LLE
x <- train
method <- "annLLE"
metric_lle <- metricML(x, fn = y, s = s, k = k, radius = radius, method = method, invert.h = TRUE, eps = 0,
                                annmethod = annmethod, distance = distance, treetype = treetype,
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

# opt.method <- c("AMISE", "MEAN", "SCALED")[2]
# riem.scale <- .1
flle <- vkde(x = fn, h = h_hdr_meta, vh = Rn, r = r, gridsize = gridsize, eval.points = fn, opt.method = opt.method, riem.scale = riem.scale, adj_matrix = adj_matrix)
## ---- echo = F------------------------------------------------------------------
p_lle <- plot_outlier(x = metric_lle, gridsize = gridsize, prob = prob, noutliers = noutliers, opt.method = opt.method, riem.scale = riem.scale, f = flle, ell.size = 0)
p_hdr_lle <- hdrscatterplot_new(E1, E2, kde.package = "ks", noutliers = noutliers, levels = prob)
# p_hdr_lle_p <- p_hdr_lle + 
#   plot_ellipse(metric_lle, ell.no = 50, add = T)
(p_lle$p + p_hdr_lle$p) + coord_fixed() + 
  plot_annotation(title = "Left: DC-KDE; Right: KDE", theme = theme(plot.title = element_text(hjust = 0.5)))
# cor(trueden, p_lle$densities, method = cormethod)
# cor(trueden, p_hdr_lle$densities, method = cormethod)

# # ellipse plot for slides
# fn <- metric_isomap$embedding
# p_isomap_ell <- fn %>%
#   as_tibble() %>%
#   ggplot(aes(x = E1, y = E2)) +
#   geom_point(color = grey(0.5), size = 1) +
#   plot_ellipse(metric_isomap, add = T, ell.no = 500, ell.size = 0,
#                color = blues9[5], fill = blues9[5], alpha = 0.2)
# fn <- metric_lle$embedding
# p_lle_ell <- fn %>%
#   as_tibble() %>%
#   ggplot(aes(x = E1, y = E2)) +
#   geom_point(color = grey(0.5), size = 1) +
#   plot_ellipse(metric_isomap, add = T, ell.no = 500, ell.size = 0,
#                color = blues9[5], fill = blues9[5], alpha = 0.2)
# p_isomap_ell + p_lle_ell






# Laplacian Eigemaps

## ----le, message=FALSE, warning=FALSE, eval=TRUE------------------------------
x <- train
method <- "annLaplacianEigenmaps"
metric_le <- metricML(x, fn = y, s = s, k = k, radius = radius, method = method, invert.h = TRUE, eps = 0,
                      annmethod = annmethod, distance = distance, treetype = treetype,
                      searchtype = searchtype
)
fn <- metric_le$embedding
Rn <- metric_le$rmetric
E1 <- fn[,1]; E2 <- fn[,2]
# fle <- vkde2d(x = E1, y = E2, h = Rn*riem.scale, gridsize = gridsize)
# fxy_le <- hdrcde:::interp.2d(fle$x, fle$y, fle$z, x0 = E1, y0 = E2)
# plot_embedding(metric_le)
# plot_ellipse(metric_le, ell.no = 50)
# plot_contour(metric_le, gridsize = gridsize, riem.scale = 1/20)

# opt.method <- c("AMISE", "MEAN", "SCALED")[3]
# riem.scale <- 1
fle <- vkde(x = fn, h = h_hdr_meta, vh = Rn, r = r, gridsize = gridsize, eval.points = fn, opt.method = opt.method, riem.scale = riem.scale, adj_matrix = adj_matrix)
## ---- echo = F------------------------------------------------------------------
p_le <- plot_outlier(x = metric_le, gridsize = gridsize, prob = prob, noutliers = noutliers, opt.method = opt.method, riem.scale = riem.scale, f = fle, ell.size = 0)
p_hdr_le <- hdrscatterplot_new(E1, E2, kde.package = "ks", noutliers = noutliers, levels = prob)
# p_hdr_le_p <- p_hdr_le + 
#   plot_ellipse(metric_le, ell.no = 50, add = T)
(p_le$p + p_hdr_le$p) + coord_fixed() + 
  plot_annotation(title = "Left: DC-KDE; Right: KDE", theme = theme(plot.title = element_text(hjust = 0.5)))
cor(trueden, p_le$densities, method = cormethod)
cor(trueden, p_hdr_le$densities, method = cormethod)


# tSNE

## ----tsne, message=FALSE, warning=FALSE, eval=TRUE------------------------------
x <- train
method <- "anntSNE"
perplexity <- 30 # round(k / 3) # 30 by default
theta <- 0 # for exact tSNE in the C++ tSNE Barnes-Hut implementation
# tictoc::tic()
metric_tsne <- metricML(x, fn = y, s = s, k = k, radius = radius, method = method, 
                        annmethod = annmethod, 
                        eps = 0, distance = distance, treetype = treetype,
                        searchtype = searchtype, 
                        perplexity = perplexity, theta = theta, invert.h = TRUE)
fn <- metric_tsne$embedding
Rn <- metric_tsne$rmetric
E1 <- fn[,1]; E2 <- fn[,2]
# ftsne <- vkde2d(x = E1, y = E2, h = Rn*riem.scale, gridsize = gridsize)
# fxy_tsne <- hdrcde:::interp.2d(ftsne$x, ftsne$y, ftsne$z, x0 = E1, y0 = E2)
# plot_embedding(metric_tsne)
# plot_ellipse(metric_tsne, ell.no = 50)
# plot_contour(metric_tsne, gridsize = gridsize, riem.scale = riem.scale)

# opt.method <- c("AMISE", "MEAN", "SCALED")[2]
# riem.scale <- .1
ftsne <- vkde(x = fn, h = h_hdr_meta, vh = Rn, r = r, gridsize = gridsize, eval.points = fn, opt.method = opt.method, riem.scale = riem.scale, adj_matrix = adj_matrix)

## ---- echo = F------------------------------------------------------------------
p_tsne <- plot_outlier(x = metric_tsne, gridsize = gridsize, prob = prob, noutliers = noutliers, opt.method = opt.method, riem.scale = riem.scale, f = ftsne, ell.size = 0)
p_hdr_tsne <- hdrscatterplot_new(E1, E2, kde.package = "ks", noutliers = noutliers, levels = prob)
# p_hdr_tsne_p <- p_hdr_tsne$p + 
#   plot_ellipse(metric_tsne, ell.no = 50, add = T)
(p_tsne$p + p_hdr_tsne$p) + coord_fixed() + 
  plot_annotation(title = "Left: DC-KDE; Right: KDE", theme = theme(plot.title = element_text(hjust = 0.5)))

# metric_tsne$embedding <- preswissroll
# plot_outlier(x = metric_tsne, gridsize = gridsize, prob = prob, noutliers = noutliers, riem.scale = 1/20, f = ftsne, ell.size = 0)
# cor(trueden, p_tsne$densities)
# cor(trueden, p_hdr_tsne$densities)


# UMAP

## ----umap, message=FALSE, warning=FALSE, eval=TRUE------------------------------
x <- train
method <- "annUMAP"
metric_umap <- metricML(x, fn = y, s = s, k = k, radius = radius, method = method, 
                        annmethod = annmethod, eps = 0, distance = distance, treetype = treetype,
                        searchtype = searchtype, 
                        invert.h = TRUE)

## ---- message=FALSE, eval=TRUE--------------------------------------------------
fn <- metric_umap$embedding
E1 <- fn[,1]; E2 <- fn[,2]
Rn <- metric_umap$rmetric

# opt.method <- c("AMISE", "MEAN", "SCALED")[2]
# riem.scale <- .1
fumap <- vkde(x = fn, h = h_hdr_meta, vh = Rn, r = r, gridsize = gridsize, eval.points = fn, opt.method = opt.method, riem.scale = riem.scale, adj_matrix = adj_matrix)
# fxy_umap <- hdrcde:::interp.2d(fumap$x, fumap$y, fumap$z, x0 = E1, y0 = E2)
# plot_embedding(metric_umap)
# plot_ellipse(metric_umap, ell.no = 50)
# plot_contour(metric_umap, gridsize = gridsize, riem.scale = 1/20)

## ---- echo = F------------------------------------------------------------------
p_umap <- plot_outlier(x = metric_umap, gridsize = gridsize, prob = prob, noutliers = noutliers, opt.method = opt.method, riem.scale = riem.scale, ell.size = 0, f = fumap)
p_hdr_umap <- hdrscatterplot_new(E1, E2, kde.package = "ks", noutliers = noutliers, levels = prob)
# p_hdr_umap_p <- p_hdr_umap +
#   plot_ellipse(metric_umap, ell.no = 50, add = T)
(p_umap$p + p_hdr_umap$p) + coord_fixed() + 
  plot_annotation(title = "Left: DC-KDE; Right: KDE", theme = theme(plot.title = element_text(hjust = 0.5)))
# cor(trueden, p_umap$densities, method = "s")
# cor(trueden, p_hdr_umap$densities, method = "s")



## -------trueoutliers------------------------------------------------------------------------
trueden <- den_2dmanifold  # NOT preswissroll$den
trueoutliers <- head(order(trueden), noutliers)
p_den_isomap <- plot_embedding_color_levels(p = plot_embedding(metric_isomap), x = trueden, N, prob, noutliers, fn = metric_isomap$embedding)
p_den_lle <- plot_embedding_color_levels(p = plot_embedding(metric_lle), x = trueden, N, prob, noutliers, fn = metric_lle$embedding)
p_den_le <- plot_embedding_color_levels(p = plot_embedding(metric_le), x = trueden, N, prob, noutliers, fn = metric_le$embedding)
p_den_tsne <- plot_embedding_color_levels(p = plot_embedding(metric_tsne), x = trueden, N, prob, noutliers, fn = metric_tsne$embedding)
p_den_umap <- plot_embedding_color_levels(p = plot_embedding(metric_umap), x = trueden, N, prob, noutliers, fn = metric_umap$embedding)

# p_den_isomap <- plot_embedding(metric_isomap) + 
#   geom_point(aes(col = trueden)) + 
#   scale_color_viridis(option = "inferno") +
#   labs(color = "Density") +
#   ggplot2::annotate("text",
#                     x = metric_isomap$embedding[outliers, 1] + diff(range(metric_isomap$embedding[,1])) / 50, 
#                     y = metric_isomap$embedding[outliers, 2] + diff(range(metric_isomap$embedding[,2])) / 50,
#                     label = outliers, col = "blue", cex = 2.5 ) 
# p_den_lle <- plot_embedding(metric_lle) + 
#   geom_point(aes(col = trueden))  + 
#   scale_color_viridis(option = "inferno") +
#   labs(color = "Density") +
#   ggplot2::annotate("text",
#                     x = metric_lle$embedding[outliers, 1] + diff(range(metric_lle$embedding[,1])) / 50, 
#                     y = metric_lle$embedding[outliers, 2] + diff(range(metric_lle$embedding[,2])) / 50,
#                     label = outliers, col = "blue", cex = 2.5 ) 
# p_den_le <- plot_embedding(metric_le) + 
#   geom_point(aes(col = trueden)) + 
#   scale_color_viridis(option = "inferno") +
#   labs(color = "Density")  +
#   ggplot2::annotate("text",
#                     x = metric_le$embedding[outliers, 1] + diff(range(metric_le$embedding[,1])) / 50, 
#                     y = metric_le$embedding[outliers, 2] + diff(range(metric_le$embedding[,2])) / 50,
#                     label = outliers, col = "blue", cex = 2.5 ) 
# p_den_tsne <- plot_embedding(metric_tsne) + 
#   geom_point(aes(col = trueden)) + 
#   scale_color_viridis(option = "inferno") +
#   labs(color = "Density") +
#   ggplot2::annotate("text",
#                     x = metric_tsne$embedding[outliers, 1] + diff(range(metric_tsne$embedding[,1])) / 50, 
#                     y = metric_tsne$embedding[outliers, 2] + diff(range(metric_tsne$embedding[,2])) / 50,
#                     label = outliers, col = "blue", cex = 2.5 ) 
# p_den_umap <- plot_embedding(metric_umap) + 
#   geom_point(aes(col = trueden)) + 
#   scale_color_viridis(option = "inferno") +
#   labs(color = "Density") +
#   ggplot2::annotate("text",
#                     x = metric_umap$embedding[outliers, 1] + diff(range(metric_umap$embedding[,1])) / 50, 
#                     y = metric_umap$embedding[outliers, 2] + diff(range(metric_umap$embedding[,2])) / 50,
#                     label = outliers, col = "blue", cex = 2.5 )  


## FINAL plot to use
# p <- (
#   (((p_den_isomap + ggtitle("ISOMAP")) | (p_den_lle + ggtitle("LLE")) | (p_den_tsne + ggtitle("t-SNE")) | (p_den_umap + ggtitle("UMAP"))) & scale_color_viridis(option = "inferno") & labs(color = "Density") & coord_fixed()) /
#     (p_isomap$p | p_lle$p | p_tsne$p | p_umap$p) /
#     (p_hdr_isomap$p | p_hdr_lle$p | p_hdr_tsne$p | p_hdr_umap$p) +  ## add _p for ellipses
#   coord_fixed() +
# #     plot_annotation(subtitle = "Top: True densities from Gaussian mixture model;
# # Middle: Outliers using variable kernel density estimate;
# # Bottom: Outliers from `hdrcde` package with KDE;") +
#   plot_layout(guides = 'collect') ) &
#   labs(x = "", y = "") &
#   theme(legend.direction = "vertical")

nolabs <- labs(x = "", y = "")
p <- (
  ( (p_den_isomap + labs(x = "", y = "ISOMAP", title = "True density") ) |
      (p_isomap$p + labs(x = "", y = "", title = "DC-KDE") ) |
      (p_hdr_isomap$p + labs(x = "", y = "", title = "KDE") ) 
  ) /
    ( (p_den_lle + labs(x = "", y = "LLE") ) | 
        (p_lle$p + nolabs) |
        (p_hdr_lle$p + nolabs) 
    ) /
    ( (p_den_le + labs(x = "", y = "Laplacian Eigenmaps") ) | 
        (p_le$p + nolabs) |
        (p_hdr_le$p + nolabs) 
    ) /
    ( (p_den_tsne + labs(x = "", y = "t-SNE") ) | 
        (p_tsne$p + nolabs) |
        (p_hdr_tsne$p + nolabs) 
    ) / 
    ( (p_den_umap + labs(x = "", y = "UMAP") ) | 
        (p_umap$p + nolabs) |
        (p_hdr_umap$p + nolabs) 
    )
) +
  plot_layout(guides = 'collect') &
  theme(legend.direction = "horizontal", legend.position = "bottom", 
        legend.box = "horizontal",
        plot.title = element_text(hjust = 0.5, face = "bold"))
p
ggsave(paste0("~/Downloads/", mapping, N, "_5levels_outliers_comparison_5ml_", opt.method, "_riem", format(if(opt.method == "SCALED") riem.scale, decimal.mark = "_"), ".png"), p, width = 8, height = 10, dpi = 300)

# ISOMAP only for slides
((p_den_isomap + labs(x = "", y = "ISOMAP", title = "True density") ) |
    (p_isomap$p + labs(x = "", y = "", title = "DC-KDE") ) |
    (p_hdr_isomap$p + labs(x = "", y = "", title = "KDE") )) +
  plot_layout(guides = 'collect') &
  theme(legend.direction = "horizontal", legend.position = "bottom", 
        legend.box = "horizontal",
        plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(paste0("~/Downloads/", mapping, N,"_5levels_outliers_comparison_isomap_3cases_", opt.method, "_riem", format(if(opt.method == "SCALED") riem.scale, decimal.mark = "_"), ".png"), width = 12, height = 5, dpi = 300)




## ----compareDensity---------------------------------------------------------------------------
fxy <- trueden
methods <- c("isomap", "lle", "le", "tsne", "umap")
# fnames <- c("fxy_vkde", "fxy_hdr")
# for(i in 1:length(methods)) {
#   phrase <- paste0("fxy_vkde <- p_", methods[1], "$densities")
#   phrase1 <- paste0("fxy_hdr <- p_hdr_", methods[1], "$densities")
#   eval(phrase)
#   eval(phrase1)
#   
# }

# # MODIFY these two lines for different ML method
# p_vkde <- p_isomap
# p_hdr <- p_hdr_isomap

## scatterplot to compare f_xy for ISOMAP
f <- tibble(fxy = fxy, # true densities
            fxy_vkde = p_isomap$densities, fxy_hdr = p_hdr_isomap$densities, # estimated densities with ISOMAP, best
            fxy_vkde_lle = p_lle$densities, fxy_hdr_lle = p_hdr_lle$densities,
            fxy_vkde_le = p_le$densities, fxy_hdr_le = p_hdr_le$densities,
            fxy_vkde_tsne = p_tsne$densities, fxy_hdr_tsne = p_hdr_tsne$densities, # estimated densities with tSNE, worst
            fxy_vkde_umap = p_umap$densities, fxy_hdr_umap = p_hdr_umap$densities
) #%>% summarise_all(rank) # Comment this part if plotting the density estimates instead of ranks
f
cor(f$fxy_vkde, f$fxy, method = "spearman")
cor(f$fxy_hdr, f$fxy, method = "spearman")
# # Plot ranks instead of densities
f <- f %>% summarise_all(rank)

pf_vkde <- f %>% 
  ggplot(aes(x = fxy, y = fxy_vkde, col = factor(preswissroll$label), shape = factor(preswissroll$label))) + 
  geom_point() + 
  labs(x = "", y = "", color = "Kernels", shape = "Kernels", title = paste("DC-KDE correlation", round(cor(f$fxy_vkde, f$fxy, method = "spearman"), 3)) ) +
  scale_y_continuous(n.breaks = 6)
pf_hdr <- f %>% 
  ggplot(aes(x = fxy, y = fxy_hdr, col = factor(preswissroll$label), shape = factor(preswissroll$label))) + 
  geom_point() + 
  labs(x = "", y = "", color = "Kernels", shape = "Kernels", title = paste("KDE correlation", round(cor(f$fxy_hdr, f$fxy, method = "spearman"), 3))) +
  scale_y_continuous(limits = c(0, max(f$fxy_hdr)), n.breaks = 5)

pf_vkde_lle <- f %>%
  ggplot(aes(x = fxy, y = fxy_vkde_lle, col = factor(preswissroll$label), shape = factor(preswissroll$label))) +
  geom_point() +
  labs(x = "", y = "", color = "Kernels", shape = "Kernels", title = paste("DC-KDE correlation", round(cor(f$fxy_vkde_lle, f$fxy, method = "spearman"), 3)) ) +
  scale_y_continuous(n.breaks = 6)
pf_hdr_lle <- f %>%
  ggplot(aes(x = fxy, y = fxy_hdr_lle, col = factor(preswissroll$label), shape = factor(preswissroll$label))) +
  geom_point() +
  labs(x = "", y = "", color = "Kernels", shape = "Kernels", title = paste("DC-KDE correlation", round(cor(f$fxy_hdr_lle, f$fxy, method = "spearman"), 3)) ) +
  scale_y_continuous(limits = c(0, max(f$fxy_hdr_lle)), n.breaks = 7)

pf_vkde_le <- f %>%
  ggplot(aes(x = fxy, y = fxy_vkde_le, col = factor(preswissroll$label), shape = factor(preswissroll$label))) +
  geom_point() +
  labs(x = "", y = "", color = "Kernels", shape = "Kernels", title = paste("DC-KDE correlation", round(cor(f$fxy_vkde_le, f$fxy, method = "spearman"), 3)) ) +
  scale_y_continuous(n.breaks = 6)
pf_hdr_le <- f %>%
  ggplot(aes(x = fxy, y = fxy_hdr_le, col = factor(preswissroll$label), shape = factor(preswissroll$label))) +
  geom_point() +
  labs(x = "", y = "", color = "Kernels", shape = "Kernels", title = paste("DC-KDE correlation", round(cor(f$fxy_hdr_le, f$fxy, method = "spearman"), 3)) ) +
  scale_y_continuous(limits = c(0, max(f$fxy_hdr_le)), n.breaks = 7)

pf_vkde_tsne <- f %>%
  ggplot(aes(x = fxy, y = fxy_vkde_tsne, col = factor(preswissroll$label), shape = factor(preswissroll$label))) +
  geom_point() +
  labs(x = "", y = "", color = "Kernels", shape = "Kernels", title = paste("DC-KDE correlation", round(cor(f$fxy_vkde_tsne, f$fxy, method = "spearman"), 3)) ) +
  scale_y_continuous(n.breaks = 6)
pf_hdr_tsne <- f %>%
  ggplot(aes(x = fxy, y = fxy_hdr_tsne, col = factor(preswissroll$label), shape = factor(preswissroll$label))) +
  geom_point() +
  labs(x = "", y = "", color = "Kernels", shape = "Kernels", title = paste("DC-KDE correlation", round(cor(f$fxy_hdr_tsne, f$fxy, method = "spearman"), 3)) ) +
  scale_y_continuous(limits = c(0, max(f$fxy_hdr_tsne)), n.breaks = 7)

pf_vkde_umap <- f %>%
  ggplot(aes(x = fxy, y = fxy_vkde_umap, col = factor(preswissroll$label), shape = factor(preswissroll$label))) +
  geom_point() +
  labs(x = "", y = "", color = "Kernels", shape = "Kernels", title = paste("DC-KDE correlation", round(cor(f$fxy_vkde_umap, f$fxy, method = "spearman"), 3)) ) +
  scale_y_continuous(n.breaks = 6)
pf_hdr_umap <- f %>%
  ggplot(aes(x = fxy, y = fxy_hdr_umap, col = factor(preswissroll$label), shape = factor(preswissroll$label))) +
  geom_point() +
  labs(x = "", y = "", color = "Kernels", shape = "Kernels", title = paste("DC-KDE correlation", round(cor(f$fxy_hdr_umap, f$fxy, method = "spearman"), 3)) ) +
  scale_y_continuous(limits = c(0, max(f$fxy_hdr_umap)), n.breaks = 7)


result <- 
  ( (pf_vkde + ylab("ISOMAP") + pf_hdr)/ 
      (pf_vkde_lle + ylab("LLE") + pf_hdr_lle) / 
      (pf_vkde_le + ylab("Laplacian Eigenmaps") + pf_hdr_le) / 
      (pf_vkde_tsne + ylab("t-SNE") + pf_hdr_tsne) / 
      (pf_vkde_umap + ylab("UMAP") + pf_hdr_umap) 
  ) +
  # (pf_vkde + pf_hdr ) + # + ylim(0, max(f[,2:3]))
  plot_layout(guides = 'collect') & 
  theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5, face = "bold"))
gt <- patchwork::patchworkGrob(result)
gt <- gridExtra::grid.arrange(gt, left = "Estimated density", bottom = "True density")

ggsave(paste0("~/Downloads/", mapping, N,"_density_comparison_5ml_riem", format(if(opt.method == "SCALED") riem.scale, decimal.mark = "_"), "_rank.png"), gt, width = 8, height = 10, dpi = 300)




# Table for density correlations
fxy <- trueden
dencor <- function(x) cor(x$densities, fxy, method = "spearman")
# dencor(p_lle)
# dencor(p_hdr_lle)
cors <- cbind(
  c(dencor(p_isomap), dencor(p_hdr_isomap)),
  c(dencor(p_lle), dencor(p_hdr_lle)),
  c(dencor(p_le), dencor(p_hdr_le)),
  c(dencor(p_tsne), dencor(p_hdr_tsne)),
  c(dencor(p_umap), dencor(p_hdr_umap))
) 
rownames(cors) <- c("DC-KDE", "KDE")
colnames(cors) <- c("ISOMAP", "LLE", "Laplacian Eigenmaps", "t-SNE", "UMAP")
cors %>% 
  kableExtra::kbl(caption = "Correlation between true density and estimated density for four manifold learning embeddings.", booktabs = TRUE, digits = 3, escape = FALSE) %>%
  kable_styling(latex_options = "scale_down") %>%
  kable_paper(full_width = TRUE)

save(cors, file = paste0("~/Downloads/CorrelationTable_", mapping, N,"_5ml_", opt.method, "_riem", format(if(opt.method == "SCALED") riem.scale, decimal.mark = "_"), ".rda"))





# Plot of outlier rank to show that the ranks do not change across different ML methods
f
prank_dckde_isolle <- f %>% 
  ggplot(aes(x = fxy_vkde, y = fxy_vkde_lle, col = factor(preswissroll$label), shape = factor(preswissroll$label))) + 
  geom_point() + 
  labs(x = "", y = "", color = "Kernels", shape = "Kernels", title = paste("DC-KDE correlation", round(cor(f$fxy_vkde, f$fxy_vkde_lle, method = "spearman"), 3)) ) +
  scale_y_continuous(n.breaks = 6)
prank_hdr_isolle <- f %>% 
  ggplot(aes(x = fxy_hdr, y = fxy_hdr_lle, col = factor(preswissroll$label), shape = factor(preswissroll$label))) + 
  geom_point() + 
  labs(x = "", y = "", color = "Kernels", shape = "Kernels", title = paste("KDE correlation", round(cor(f$fxy_hdr, f$fxy_hdr_lle, method = "spearman"), 3))) +
  scale_y_continuous(limits = c(0, max(f$fxy_hdr)), n.breaks = 5)

prank_dckde_isole <- f %>% 
  ggplot(aes(x = fxy_vkde, y = fxy_vkde_le, col = factor(preswissroll$label), shape = factor(preswissroll$label))) + 
  geom_point() + 
  labs(x = "", y = "", color = "Kernels", shape = "Kernels", title = paste("DC-KDE correlation", round(cor(f$fxy_vkde, f$fxy_vkde_le, method = "spearman"), 3)) ) +
  scale_y_continuous(n.breaks = 6)
prank_hdr_isole <- f %>% 
  ggplot(aes(x = fxy_hdr, y = fxy_hdr_le, col = factor(preswissroll$label), shape = factor(preswissroll$label))) + 
  geom_point() + 
  labs(x = "", y = "", color = "Kernels", shape = "Kernels", title = paste("KDE correlation", round(cor(f$fxy_hdr, f$fxy_hdr_le, method = "spearman"), 3))) # +
  # scale_y_continuous(limits = c(0, max(f$fxy_hdr)), n.breaks = 5)

prank_dckde_isotsne <- f %>% 
  ggplot(aes(x = fxy_vkde, y = fxy_vkde_tsne, col = factor(preswissroll$label), shape = factor(preswissroll$label))) + 
  geom_point() + 
  labs(x = "", y = "", color = "Kernels", shape = "Kernels", title = paste("DC-KDE correlation", round(cor(f$fxy_vkde, f$fxy_vkde_tsne, method = "spearman"), 3)) ) +
  scale_y_continuous(n.breaks = 6)
prank_hdr_isotsne <- f %>% 
  ggplot(aes(x = fxy_hdr, y = fxy_hdr_tsne, col = factor(preswissroll$label), shape = factor(preswissroll$label))) + 
  geom_point() + 
  labs(x = "", y = "", color = "Kernels", shape = "Kernels", title = paste("KDE correlation", round(cor(f$fxy_hdr, f$fxy_hdr_tsne, method = "spearman"), 3))) # +
  # scale_y_continuous(limits = c(0, max(f$fxy_hdr)), n.breaks = 5)

prank_dckde_isoumap <- f %>% 
  ggplot(aes(x = fxy_vkde, y = fxy_vkde_umap, col = factor(preswissroll$label), shape = factor(preswissroll$label))) + 
  geom_point() + 
  labs(x = "", y = "", color = "Kernels", shape = "Kernels", title = paste("DC-KDE correlation", round(cor(f$fxy_vkde, f$fxy_vkde_umap, method = "spearman"), 3)) ) +
  scale_y_continuous(n.breaks = 6)
prank_hdr_isoumap <- f %>% 
  ggplot(aes(x = fxy_hdr, y = fxy_hdr_umap, col = factor(preswissroll$label), shape = factor(preswissroll$label))) + 
  geom_point() + 
  labs(x = "", y = "", color = "Kernels", shape = "Kernels", title = paste("KDE correlation", round(cor(f$fxy_hdr, f$fxy_hdr_umap, method = "spearman"), 3))) # +
  # scale_y_continuous(limits = c(0, max(f$fxy_hdr)), n.breaks = 5)

(prank_dckde_isolle + prank_hdr_isolle)

prank_isomap <- (
    (prank_dckde_isolle + ylab("LLE") + prank_hdr_isolle) / 
    (prank_dckde_isole + ylab("Laplacian Eigenmaps") + prank_hdr_isole) / 
    (prank_dckde_isotsne + ylab("t-SNE") + prank_hdr_isotsne) / 
    (prank_dckde_isoumap + ylab("UMAP") + prank_hdr_isoumap) 
) +
  # (pf_vkde + pf_hdr ) + # + ylim(0, max(f[,2:3]))
  plot_layout(guides = 'collect') & 
  theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5, face = "bold"))
prank_isomap
ggsave(paste0("~/Downloads/", mapping, N,"_density_comparison_isomapvs4ml_riem", format(if(opt.method == "SCALED") riem.scale, decimal.mark = "_"), "_rank.png"), prank_isomap, width = 8, height = 10, dpi = 300)




# TODO:
# Spec and sensitivity for correctly categorized outliers, over a small area in the top left
# Change plot of true density colors to 5 levels
# Table of 
# 
