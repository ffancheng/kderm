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
N <- 2000
p <- 2

##----Data Generation-------------------------------------
## 1. Generate 2-D manifold and 3-D twin peaks
##--------------------------------------------------------
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



## -------------------------------------------------------------------------------
# 2. TRUE density of 2-D manifold
## -------------------------------------------------------------------------------
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
                          annmethod = annmethod, distance = distance, treetype = treetype,
                          searchtype = searchtype
)
# summary(metric_meta)
# all.equal(y, metric_meta$embedding, check.attributes = F)
fn <- metric_meta$embedding
Rn <- metric_meta$rmetric

# ## ----ggellipse, include=FALSE, eval=FALSE-----------------------------------
# plot_embedding(metric_meta) +
#   labs(x = "ISO1", y = "ISO2")
# plot_ellipse(metric_meta, add = F, ell.no = 500, ell.size = .1,
#              color = blues9[5], fill = blues9[5], alpha = 0.2)


# Run once with y <- sr$metadata
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
# ggsave("paper/figures/truedensity_twinpeaks_dc.png", plotmanifold, width = 8, height = 6, dpi = 300)
plot_ly(data = swissroll, x = ~ x, y = ~ y, z = ~ z, color = den_2dmanifold, # colored with 3d density
        type = "scatter3d", mode = "markers", size = 1, text = paste("density:", preswissroll$den))
trueden <- den_2dmanifold # sr$den
trueoutliers <- head(order(trueden), noutliers)
# save(sr, train, trueden, den_twinpeaks, noutliers, trueoutliers, file = paste0("data/simdata_3d_N", N, "_01_trueden_k", k, ".rda"))



## -------------------------------------------------------------------------------
# 3. THE FOLLOWING ESTIMATION IS BASED ON THE META DATA (Not run. To check if dc-kde works)
## -------------------------------------------------------------------------------
# KDE with KDE
E1 <- fn[,1] # rename as Ed to match the aesthetics in plot_ellipse()
E2 <- fn[,2]
prob <- c(1, 10, 50, 95, 99) #
p_hdr_meta <- hdrscatterplot_new(E1, E2, kde.package = "ks", levels = prob, noutliers = noutliers, label = NULL)
h_hdr_meta <- p_hdr_meta$den$den$h # same as ks::Hpi.diag(fn)
p_hdr_meta$densities %>% summary()

# DC-KDE
adj_matrix <- metric_meta$adj_matrix
opt.method <- c("AMISE", "MEAN", "SCALED")[3] # 3 ONLY FOR NOW, no scaling for Rn
riem.scale <- 1 # h_hdr_meta
r <- .5 # 0.829 vs 0.773
# r <- sqrt(mean(apply(Rn, 3, det)))
fmeta <- vkde(x = fn, h = h_hdr_meta, vh = Rn, r = r, gridsize = gridsize, eval.points = fn, opt.method = opt.method, riem.scale = riem.scale, adj_matrix = adj_matrix)
p_meta <- plot_outlier(x = metric_meta, gridsize = gridsize, prob = prob, riem.scale = riem.scale, f = fmeta, ell.size = 0)
p_meta$densities %>% summary()

## TODO: save computed results, beautify plots, remove same y axis, scale_x_log to show low density area, edit legend to one line, change color in rank plot to categorized HDR colors


# Rank correlation
f <- tibble(fxy = trueden, # true densities
            fxy_dckde = p_meta$densities, 
            fxy_kde = p_hdr_meta$densities
)
cormethod <- c("pearson", "kendall", "spearman")[3]
(rankcompare_meta <- 
c("dckde cor" = cor(f$fxy_dckde, f$fxy, method = cormethod), 
  "kde cor" = cor(f$fxy_kde, f$fxy, method = cormethod),
  "dckde MSE" = mean((f$fxy_dckde - f$fxy) ^ 2),
  "kde MSE" = mean((f$fxy_kde - f$fxy) ^ 2)
))

# # Plot ranks instead of densities
frank <- f %>% summarise_all(rank)

# save(metric_meta, p_hdr_meta, r, fmeta, p_meta, f, frank, rankcompare_meta, file = paste0("data/simdata_3d_N", N, "_02_truemanifold_k", k, "_dckde_r", format(r, decimal.mark = "_"), ".rda"))


## -----------------------------------------------------------------------------
# # Plotting (Not run)
# pf_vkde <- frank %>% 
#   ggplot(aes(x = fxy, y = fxy_dckde, col = factor(preswissroll$label), shape = factor(preswissroll$label))) + 
#   geom_point() + 
#   labs(x = "", y = "", color = "Kernels", shape = "Kernels", title = paste("DC-KDE correlation", format(cor(f$fxy_dckde, f$fxy, method = "spearman"), 3)) ) +
#   scale_y_continuous(n.breaks = 5)
# pf_hdr <- frank %>% 
#   ggplot(aes(x = fxy, y = fxy_kde, col = factor(preswissroll$label), shape = factor(preswissroll$label))) + 
#   geom_point() + 
#   labs(x = "", y = "", color = "Kernels", shape = "Kernels", title = paste("KDE correlation", format(cor(f$fxy_kde, f$fxy, method = "spearman"), 3))) +
#   scale_y_continuous(n.breaks = 5)
# # Rank plot
# pf_vkde + pf_hdr

trueoutliers <- head(order(trueden), noutliers)
frank_long <- frank %>% 
  mutate(outlier = case_when((fxy <= noutliers) ~ "True outliers",
                             (fxy_dckde <= noutliers) ~ "DC-KDE outliers",
                             (fxy_kde <= noutliers) ~ "KDE outliers",
                             TRUE ~ "Not outliers")) %>% 
  pivot_longer(cols = -c(fxy, outlier), names_to = "kde", values_to = "densities")
frank_long %>% 
  ggplot(aes(fxy, densities, col = outlier)) + 
  geom_point() +
  facet_wrap(~ kde, 
             labeller = as_labeller(c(
               `fxy_dckde` = paste("DC-KDE correlation:", format(rankcompare_meta[1], digits = 3)),
               `fxy_kde` = paste("KDE correlation:", format(rankcompare_meta[2], digits = 3))))
             ) + 
  scale_x_log10() + 
  scale_y_log10() +
  labs(x = "True density rank", y = "Estimated density rank")
  # annotate("text",
  #          x = frank_long[frank_long$fxy %in% trueoutliers, "fxy"] + 1, 
  #          y = frank_long[frank_long$fxy %in% trueoutliers, "densities"] + 1,
  #          label = paste(trueoutliers), col = "blue", cex = 2.5)


# p_den_meta <- plot_embedding(metric_meta) +
#   geom_point(aes(col = trueden)) +
#   scale_color_viridis(option = "inferno") +
#   labs(color = "Density") +
#   ggplot2::annotate("text",
#                     x = metric_meta$embedding[trueoutliers, 1] + diff(range(metric_meta$embedding[,1])) / 50,
#                     y = metric_meta$embedding[trueoutliers, 2] + diff(range(metric_meta$embedding[,2])) / 50,
#                     label = trueoutliers, col = "blue", cex = 2.5 )

### Categorize density to 5 levels prob
# plot_embedding_color_levels <- function(x = trueden, N, prob = c(1, 50, 95, 99), noutliers, fn = NULL) {
#   if(max(prob) > 50){ # Assume prob is coverage percentage
#     alpha <- (100 - prob) / 100
#   } else {# prob is tail probability (for backward compatibility)
#     alpha <- prob
#   }
#   alpha <- sort(alpha) # 0.99 0.50 0.05 0.01
#   falpha <- quantile(x, alpha) # increasing
#   region <- numeric(N) + 100 # rep(100, N)
#   for (i in seq_along(prob)) region[x > falpha[i]] <- 100 - prob[i]
#   if (is.null(noutliers)) noutliers <- sum(region > max(prob))
#   noutliers <- min(noutliers, N)
#   
#   # Construct region factor
#   prob <- sort(unique(region[region < 100]), decreasing = TRUE)
#   prob <- c(prob, 100)
#   Region <- factor(region,
#                    levels = prob,
#                    labels = c(paste(head(prob, -1)), ">99") # "99"  "50"  "5"   "1"   ">99"
#   )
#   # Sort data so the larger regions go first (other than outliers)
#   # rr <- region
#   # rr[region == 100] <- 0
#   # ord <- order(rr, decreasing = TRUE)
#   ord <- order(x, decreasing = T)  # outliers at the end
#   # fn_dec <- fn[ord,] # data is reordered by ord
#   outliers <- tail(ord, noutliers)
#   
#   # if (noutliers > 0) {
#   #   outlier_rank <- order(x[ord])  # order density as well to match the new data
#   #   outliers <- outlier_rank[seq_len(noutliers)] # take top noutliers labels for annotation
#   # }
# 
#   p <- fn %>% 
#     as_tibble() %>% 
#     ggplot(aes(x = E1, y = E2)) + 
#     ggplot2::geom_point(ggplot2::aes_string(col = "Region")) + 
#     ggplot2::scale_color_manual(
#       name = "HDRs",
#       breaks = c(paste(head(sort(prob), -1)), ">99"),
#       values = c(RColorBrewer::brewer.pal(length(prob), "YlOrRd")[-1], "#000000")) + 
#     ggplot2::annotate("text",
#                       x = fn[outliers, 1] + diff(range(fn[outliers, 1])) / 50, 
#                       y = fn[outliers, 2] + diff(range(fn[outliers, 2])) / 50,
#                       label = outliers, # seq_len(N)[outliers], 
#                       col = "blue", cex = 2.5
#                       )
#   return(p)
# }

p_den_meta <- plot_embedding_color_levels(x = trueden, N, prob, noutliers, fn = metric_meta$embedding)

(p_den_meta + p_meta$p + p_hdr_meta$p) + 
  coord_fixed() + 
  plot_annotation(title = "Left: True manifold density; Middle: DC-KDE; Right: KDE", 
                  theme = theme(plot.title = element_text(hjust = 0.5)))
# ggsave("~/Downloads/hdrplot_meta_metafn_density_compare.png", width = 12, height = 6)
# ggsave(paste0("~/Downloads/hdrplot_swissroll_meta_metafn_density_compare_thetaform2_riem", r, "_", format(if(opt.method == "SCALED") riem.scale, decimal.mark = "_"), ".png"), width = 12, height = 6)





## -------------------------------------------------------------------------------
# 4. THE FOLLOWING ESTIMATION IS BASED ON THE MANIFOLD LEARNING EMBEDDING (Run from here)
## -------------------------------------------------------------------------------
# 4.1 ISOMAP
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

p_den_isomap <- plot_embedding_color_levels(x = trueden, N, prob, noutliers, fn = metric_isomap$embedding)
(p_den_isomap + p_isomap$p + p_hdr_isomap$p) + 
  coord_fixed()


# 4.2 LLE
x <- train
method <- "annLLE"
metric_lle <- metricML(x, fn = y, s = s, k = k, radius = radius, method = method, invert.h = TRUE, eps = 0,
                                annmethod = annmethod, distance = distance, treetype = treetype,
                                searchtype = searchtype
)
fn <- metric_lle$embedding
Rn <- metric_lle$rmetric
adj_matrix <- metric_lle$adj_matrix
E1 <- fn[,1]; E2 <- fn[,2]
# flle <- vkde2d(x = E1, y = E2, h = Rn*riem.scale, gridsize = gridsize)
# fxy_lle <- hdrcde:::interp.2d(flle$x, flle$y, flle$z, x0 = E1, y0 = E2)
# plot_embedding(metric_lle)
# plot_ellipse(metric_lle, ell.no = 50)
# plot_contour(metric_lle, gridsize = gridsize, riem.scale = 1/20)

# opt.method <- c("AMISE", "MEAN", "SCALED")[2]
# riem.scale <- .1
flle <- vkde(x = fn, h = h_hdr_meta, vh = Rn, r = r, gridsize = gridsize, eval.points = fn, opt.method = opt.method, riem.scale = riem.scale, adj_matrix = adj_matrix)
## ----------------------------------------------------------------------
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



# 4.3 Laplacian Eigemaps
x <- train
method <- "annLaplacianEigenmaps"
metric_le <- metricML(x, fn = y, s = s, k = k, radius = radius, method = method, invert.h = TRUE, eps = 0,
                      annmethod = annmethod, distance = distance, treetype = treetype,
                      searchtype = searchtype
)
fn <- metric_le$embedding
Rn <- metric_le$rmetric
adj_matrix <- metric_le$adj_matrix
E1 <- fn[,1]; E2 <- fn[,2]
# fle <- vkde2d(x = E1, y = E2, h = Rn*riem.scale, gridsize = gridsize)
# fxy_le <- hdrcde:::interp.2d(fle$x, fle$y, fle$z, x0 = E1, y0 = E2)
# plot_embedding(metric_le)
# plot_ellipse(metric_le, ell.no = 50)
# plot_contour(metric_le, gridsize = gridsize, riem.scale = 1/20)

# opt.method <- c("AMISE", "MEAN", "SCALED")[3]
# riem.scale <- 1
fle <- vkde(x = fn, h = h_hdr_meta, vh = Rn, r = r, gridsize = gridsize, eval.points = fn, opt.method = opt.method, riem.scale = riem.scale, adj_matrix = adj_matrix)
## ---------------------------------------------------------------------
p_le <- plot_outlier(x = metric_le, gridsize = gridsize, prob = prob, noutliers = noutliers, opt.method = opt.method, riem.scale = riem.scale, f = fle, ell.size = 0)
p_hdr_le <- hdrscatterplot_new(E1, E2, kde.package = "ks", noutliers = noutliers, levels = prob)
# p_hdr_le_p <- p_hdr_le + 
#   plot_ellipse(metric_le, ell.no = 50, add = T)
(p_le$p + p_hdr_le$p) + coord_fixed() + 
  plot_annotation(title = "Left: DC-KDE; Right: KDE", theme = theme(plot.title = element_text(hjust = 0.5)))
cor(trueden, p_le$densities, method = cormethod)
cor(trueden, p_hdr_le$densities, method = cormethod)


# 4.4 tSNE
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
adj_matrix <- metric_tsne$adj_matrix
E1 <- fn[,1]; E2 <- fn[,2]
# ftsne <- vkde2d(x = E1, y = E2, h = Rn*riem.scale, gridsize = gridsize)
# fxy_tsne <- hdrcde:::interp.2d(ftsne$x, ftsne$y, ftsne$z, x0 = E1, y0 = E2)
# plot_embedding(metric_tsne)
# plot_ellipse(metric_tsne, ell.no = 50)
# plot_contour(metric_tsne, gridsize = gridsize, riem.scale = riem.scale)

# opt.method <- c("AMISE", "MEAN", "SCALED")[2]
# riem.scale <- .1
ftsne <- vkde(x = fn, h = h_hdr_meta, vh = Rn, r = r, gridsize = gridsize, eval.points = fn, opt.method = opt.method, riem.scale = riem.scale, adj_matrix = adj_matrix)

## dc-kde
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


# 4.5 UMAP
x <- train
method <- "annUMAP"
metric_umap <- metricML(x, fn = y, s = s, k = k, radius = radius, method = method, 
                        annmethod = annmethod, eps = 0, distance = distance, treetype = treetype,
                        searchtype = searchtype, 
                        invert.h = TRUE)
fn <- metric_umap$embedding
E1 <- fn[,1]; E2 <- fn[,2]
Rn <- metric_umap$rmetric
adj_matrix <- metric_umap$adj_matrix

fumap <- vkde(x = fn, h = h_hdr_meta, vh = Rn, r = r, gridsize = gridsize, eval.points = fn, opt.method = opt.method, riem.scale = riem.scale, adj_matrix = adj_matrix)
# fxy_umap <- hdrcde:::interp.2d(fumap$x, fumap$y, fumap$z, x0 = E1, y0 = E2)
# plot_embedding(metric_umap)
# plot_ellipse(metric_umap, ell.no = 50)
# plot_contour(metric_umap, gridsize = gridsize, riem.scale = 1/20)

## ----------------------------------------------------------------------
p_umap <- plot_outlier(x = metric_umap, gridsize = gridsize, prob = prob, noutliers = noutliers, opt.method = opt.method, riem.scale = riem.scale, ell.size = 0, f = fumap)
p_hdr_umap <- hdrscatterplot_new(E1, E2, kde.package = "ks", noutliers = noutliers, levels = prob)
# p_hdr_umap_p <- p_hdr_umap +
#   plot_ellipse(metric_umap, ell.no = 50, add = T)
(p_umap$p + p_hdr_umap$p) + coord_fixed() + 
  plot_annotation(title = "Left: DC-KDE; Right: KDE", theme = theme(plot.title = element_text(hjust = 0.5)))
# cor(trueden, p_umap$densities, method = "s")
# cor(trueden, p_hdr_umap$densities, method = "s")


# save(metric_isomap, fisomap, p_isomap, p_hdr_isomap,
#      metric_lle, flle, p_lle, p_hdr_lle,
#      metric_le, fle, p_le, p_hdr_le,
#      metric_tsne, ftsne, p_tsne, p_hdr_tsne,
#      metric_umap, fumap, p_umap, p_hdr_umap,
#      prob, r, radius,
#      file = paste0("data/simdata_3d_N", N, "_03_5mlembedding_radius", radius, "_dckde_r", format(r, decimal.mark = "_"), "_rda")
#      )


## ---- Plotting----------------------------------------------------------------
## 5. Embedding plots colored with density levels 
## -----------------------------------------------------------------------------
# trueden <- den_2dmanifold  # NOT preswissroll$den
# trueoutliers <- head(order(trueden), noutliers)
p_den_isomap <- plot_embedding_color_levels(x = trueden, N, prob, noutliers, fn = metric_isomap$embedding)
p_den_lle <- plot_embedding_color_levels(x = trueden, N, prob, noutliers, fn = metric_lle$embedding)
p_den_le <- plot_embedding_color_levels(x = trueden, N, prob, noutliers, fn = metric_le$embedding)
p_den_tsne <- plot_embedding_color_levels(x = trueden, N, prob, noutliers, fn = metric_tsne$embedding)
p_den_umap <- plot_embedding_color_levels(x = trueden, N, prob, noutliers, fn = metric_umap$embedding)

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
noytext <-  theme(axis.text.y = element_blank())
p <- (
  (( (p_den_isomap + labs(x = "", y = "ISOMAP", title = "True density") ) |
      (p_isomap$p + labs(x = "", y = "", title = "DC-KDE") + noytext) |
      (p_hdr_isomap$p + labs(x = "", y = "", title = "KDE") + noytext ) 
  ) & scale_y_continuous(n.breaks = 4)) /
    (( (p_den_lle + labs(x = "", y = "LLE") ) | 
        (p_lle$p + nolabs + noytext) |
        (p_hdr_lle$p + nolabs + noytext) 
    ) & scale_y_continuous(breaks = c(-5, 0, 5))) /
    (( (p_den_le + labs(x = "", y = "Laplacian Eigenmaps") ) | 
        (p_le$p + nolabs + noytext) |
        (p_hdr_le$p + nolabs + noytext) 
    ) & scale_x_continuous(breaks = c(-.025, 0, .025)) & scale_y_continuous(breaks = c(-.05, 0, .05), limits = c(-.05, .05)) ) /
    (( (p_den_tsne + labs(x = "", y = "t-SNE") ) | 
        (p_tsne$p + nolabs + noytext) |
        (p_hdr_tsne$p + nolabs + noytext) 
    ) & scale_y_continuous(n.breaks = 4)) / 
    (( (p_den_umap + labs(x = "", y = "UMAP") ) | 
        (p_umap$p + nolabs + noytext) |
        (p_hdr_umap$p + nolabs + noytext) 
    ) & scale_y_continuous(n.breaks = 4))
) + 
  plot_layout(guides = 'collect') &
  # guides(HDRs = guide_legend(nrow = 1)) &
  # scale_y_continuous(n.breaks = 4) &
  theme(legend.direction = "horizontal", legend.position = "bottom", 
        legend.box = "horizontal", axis.title.y = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"))
p
# `tpoutliers` chunk
# ggsave(paste0("paper/figures/", mapping, N, "_5levels_outliers_comparison_5ml_r", format(r, decimal.mark = "_"), ".png"), p, width = 8, height = 10, dpi = 300)


# ISOMAP only for slides
((p_den_isomap + labs(x = "", y = "ISOMAP", title = "True density") ) |
    (p_isomap$p + labs(x = "", y = "", title = "DC-KDE") + noytext) |
    (p_hdr_isomap$p + labs(x = "", y = "", title = "KDE") + noytext)) +
  plot_layout(guides = 'collect') &
  scale_x_continuous(n.breaks = 4) &
  scale_y_continuous(n.breaks = 4) &
  theme(legend.direction = "horizontal", legend.position = "bottom", 
        legend.box = "horizontal",
        plot.title = element_text(hjust = 0.5, face = "bold"))
# ggsave(paste0("paper/figures/", mapping, N,"_5levels_outliers_comparison_isomap_3cases_r", format(r, decimal.mark = "_"), ".png"), width = 12, height = 5, dpi = 300)



## -----------------------------------------------------------------------------
## 6. Plots and tables to compare density ranks
## -----------------------------------------------------------------------------
fxy <- trueden
methods <- c("isomap", "lle", "le", "tsne", "umap")

## scatterplot to compare f_xy for ISOMAP
f <- tibble(fxy = fxy, # true densities
            fxy_vkde = p_isomap$densities, fxy_hdr = p_hdr_isomap$densities, # estimated densities with ISOMAP, best
            fxy_vkde_lle = p_lle$densities, fxy_hdr_lle = p_hdr_lle$densities,
            fxy_vkde_le = p_le$densities, fxy_hdr_le = p_hdr_le$densities,
            fxy_vkde_tsne = p_tsne$densities, fxy_hdr_tsne = p_hdr_tsne$densities, # estimated densities with tSNE, worst
            fxy_vkde_umap = p_umap$densities, fxy_hdr_umap = p_hdr_umap$densities
)
f
# # Plot ranks instead of densities
frank <- f %>% summarise_all(rank)

## -----------------------------------------------------------------------------
# 6.1 Table for density correlations
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
  kable_paper(full_width = TRUE) %>%
  row_spec(1, bold = TRUE)

# save(cors, file = paste0("paper/figures/CorrelationTable_", mapping, N,"_5ml_radius", radius, "_r", format(r, decimal.mark = "_"), ".rda"))



## -----------------------------------------------------------------------------
# 6.2 Rank comparison plots, dckde vs kde
frank_long <- frank %>%
  # mutate(outlier = case_when((fxy <= noutliers) ~ "True outliers",
  #                            (fxy_dckde <= noutliers) ~ "DC-KDE outliers",
  #                            (fxy_kde <= noutliers) ~ "KDE outliers",
  #                            TRUE ~ "Not outliers")) %>%
  rename(fxy_vkde_isomap = fxy_vkde,
         fxy_hdr_isomap = fxy_hdr) %>%
  `colnames<-`(gsub("fxy_", "", colnames(.))) %>% 
  pivot_longer(cols = -c(fxy), names_to = "kde", values_to = "densities") %>% 
  separate(kde, c("KDE", "ML"), sep = "_") %>% 
  mutate(KDE = factor(KDE, levels = c("vkde", "hdr"), labels = c("DC-KDE", "KDE")),
         ML = factor(ML, levels = methods, labels = c("ISOMAP", "LLE", "Laplacian Eigenmaps", "t-SNE", "UMAP"))
         )
gt <- frank_long %>%
  ggplot(aes(fxy, densities)) +
  geom_point() +
  facet_grid(ML ~ KDE) +
  # scale_x_log10() +
  # scale_y_log10() +
  labs(x = "True density rank", y = "Estimated density rank")
gt
# ggsave(paste0("paper/figures/", mapping, N,"_density_comparison_5ml_radius", radius, "_r", format(r, decimal.mark = "_"), "_rank_nocor.png"), gt, width = 8, height = 10, dpi = 300)


# Not run
# pf_vkde <- frank %>% 
#   ggplot(aes(x = fxy, y = fxy_vkde)) +  # , col = factor(preswissroll$label)
#   geom_point() + 
#   labs(x = "", y = "", color = "Kernels", subtitle = paste("Rank correlation =", format(cor(f$fxy_vkde, f$fxy, method = "spearman"), digits = 3)) ) +
#   scale_y_continuous(n.breaks = 6)
# pf_hdr <- frank %>% 
#   ggplot(aes(x = fxy, y = fxy_hdr)) + 
#   geom_point() + 
#   labs(x = "", y = "", color = "Kernels", subtitle = paste("Rank correlation =", format(cor(f$fxy_hdr, f$fxy, method = "spearman"), digits = 3))) +
#   scale_y_continuous(n.breaks = 5)
# pf_vkde_lle <- frank %>%
#   ggplot(aes(x = fxy, y = fxy_vkde_lle)) +
#   geom_point() +
#   labs(x = "", y = "", color = "Kernels", subtitle = paste("Rank correlation =", format(cor(f$fxy_vkde_lle, f$fxy, method = "spearman"), digits = 3)) ) +
#   scale_y_continuous(n.breaks = 6)
# pf_hdr_lle <- frank %>%
#   ggplot(aes(x = fxy, y = fxy_hdr_lle)) +
#   geom_point() +
#   labs(x = "", y = "", color = "Kernels", subtitle = paste("Rank correlation =", format(cor(f$fxy_hdr_lle, f$fxy, method = "spearman"), digits = 3, nsmall = 3)) ) +
#   scale_y_continuous(n.breaks = 7)
# pf_vkde_le <- frank %>%
#   ggplot(aes(x = fxy, y = fxy_vkde_le)) +
#   geom_point() +
#   labs(x = "", y = "", color = "Kernels", subtitle = paste("Rank correlation =", format(cor(f$fxy_vkde_le, f$fxy, method = "spearman"), digits = 3)) ) +
#   scale_y_continuous(n.breaks = 6)
# pf_hdr_le <- frank %>%
#   ggplot(aes(x = fxy, y = fxy_hdr_le)) +
#   geom_point() +
#   labs(x = "", y = "", color = "Kernels", subtitle = paste("Rank correlation =", format(cor(f$fxy_hdr_le, f$fxy, method = "spearman"), digits = 3)) ) +
#   scale_y_continuous(n.breaks = 7)
# pf_vkde_tsne <- frank %>%
#   ggplot(aes(x = fxy, y = fxy_vkde_tsne)) +
#   geom_point() +
#   labs(x = "", y = "", color = "Kernels", subtitle = paste("Rank correlation =", format(cor(f$fxy_vkde_tsne, f$fxy, method = "spearman"), digits = 3)) ) +
#   scale_y_continuous(n.breaks = 6)
# pf_hdr_tsne <- frank %>%
#   ggplot(aes(x = fxy, y = fxy_hdr_tsne)) +
#   geom_point() +
#   labs(x = "", y = "", color = "Kernels", subtitle = paste("Rank correlation =", format(cor(f$fxy_hdr_tsne, f$fxy, method = "spearman"), digits = 3)) ) +
#   scale_y_continuous(n.breaks = 7)
# pf_vkde_umap <- frank %>%
#   ggplot(aes(x = fxy, y = fxy_vkde_umap)) +
#   geom_point() +
#   labs(x = "", y = "", color = "Kernels", subtitle = paste("Rank correlation =", format(cor(f$fxy_vkde_umap, f$fxy, method = "spearman"), digits = 3)) ) +
#   scale_y_continuous(n.breaks = 6)
# pf_hdr_umap <- frank %>%
#   ggplot(aes(x = fxy, y = fxy_hdr_umap)) +
#   geom_point() +
#   labs(x = "", y = "", color = "Kernels", subtitle = paste("Rank correlation =", format(cor(f$fxy_hdr_umap, f$fxy, method = "spearman"), digits = 3)) ) +
#   scale_y_continuous(n.breaks = 7)
#
# result <- 
#   ( (pf_vkde + labs(y = "ISOMAP", title = "DC-KDE") + (pf_hdr + labs(title = "KDE") + theme(axis.text.y = element_blank())))/ 
#       (pf_vkde_lle + labs(y = "LLE") + (pf_hdr_lle + theme(axis.text.y = element_blank()))) / 
#       (pf_vkde_le + ylab("Laplacian Eigenmaps") + (pf_hdr_le + theme(axis.text.y = element_blank()))) / 
#       (pf_vkde_tsne + ylab("t-SNE") + (pf_hdr_tsne + theme(axis.text.y = element_blank()))) / 
#       (pf_vkde_umap + ylab("UMAP") + (pf_hdr_umap + theme(axis.text.y = element_blank()))) 
#   ) +
#   plot_layout(guides = 'collect') & 
#   theme(plot.title = element_text(hjust = 0.5, face = "bold")) &  # legend.position = 'bottom'
#   geom_point(color = "black", size = .5) & 
#   guides(color = "none")
# gt <- patchwork::patchworkGrob(result)
# gt <- gridExtra::grid.arrange(gt, left = "Estimated density", bottom = "True density")
# ggsave(paste0("paper/figures/", mapping, N,"_density_comparison_5ml_radius", radius, "_r", format(r, decimal.mark = "_"), "_rank.png"), gt, width = 8, height = 10, dpi = 300)






## -----------------------------------------------------------------------------
## 6.3 Only isomap for the paper (Figure `tpisomapden`)
# RUN this!!!
##---------------------------------------------------------------------
# color the points with HDR, add rectangle and highlight points outside
xname <- c("isomap", "tsne")[2] # !!! Modify index for different ML methods comparison
probs <- prob[1:2]
levels <- c(probs / 100, 1)
q_rank <- quantile(1:N, levels)
rec <- log10(q_rank)
greys <- scales::grey_pal()(10)[c(6, 8, 10)]
frank1 <- frank %>% 
  # select(fxy:fxy_hdr) %>% # isomap
  rename(fxy_vkde_isomap = fxy_vkde,
         
         fxy_hdr_isomap = fxy_hdr) 
frank1 <- frank1 %>%
  select(fxy, contains(paste0("_", xname))) %>% 
  mutate(across(everything(), 
                function(x) cut(x, 
                                breaks = c(0, q_rank), 
                                labels = c(">99", paste(100 - probs))))
  ) %>% 
  `colnames<-`(gsub("fxy", "col", colnames(.))) %>% 
  bind_cols(select(frank1, fxy, contains(paste0("_", xname)))) %>% 
  `colnames<-`(gsub(paste0("_", xname), "", colnames(.))) %>% 
  mutate(dckde_false = (col == col_vkde) & (col != last(probs)) & (col_vkde != last(probs)),
         kde_false = (col == col_hdr) & (col != last(probs)) & (col_hdr != last(probs)) ) 

plogrec_isomap_vkde <-
  frank1 %>%
  ggplot(aes(x = fxy, y = fxy_vkde)) +
  # geom_rect(xmin = 0, xmax = rec[2], ymin = 0, ymax = rec[2], fill = greys[2], col = greys[2]) +
  # geom_rect(xmin = rec[2], xmax = rec[3], ymin = rec[2], ymax = rec[3], fill = greys[3], col = greys[3]) + # third blocks
  geom_rect(xmin = rec[1], xmax = rec[2], ymin = rec[1], ymax = rec[2], fill = greys[2], col = greys[2]) + # small blocks
  geom_rect(xmin = 0, xmax = rec[1], ymin = 0, ymax = rec[1], fill = greys[1], col = greys[1]) +
  geom_point(aes(col = col_vkde, shape = dckde_false, size = dckde_false, alpha = dckde_false)) + # if color HDRs
  scale_color_manual(
    name = "HDRs",
    breaks = levels(frank1$col),
    values = rev(c(RColorBrewer::brewer.pal(length(prob)+1, "YlOrRd")[-1], "#000000"))[1:length(levels)],
    labels = paste0(levels(frank1$col), "%")
  ) +
  scale_x_continuous(trans = "log10", limits = c(1, 2000), breaks = c(0, 20, 200, 1000, 2000), guide = guide_axis(n.dodge = 1)) + 
  scale_y_continuous(trans = "log10", limits = c(1, 2000), breaks = c(0, 20, 200, 1000, 2000)) + 
  coord_fixed() +
  scale_shape_manual(values = c(`TRUE` = 20, `FALSE` = 17 )) +
  scale_size_manual(values = c(`TRUE` = 1, `FALSE` = 2 )) +
  scale_alpha_manual(values = c(`TRUE` = .5, `FALSE` = 1) ) + 
  labs(x = "", y = "", shape = "Correct HDRs", size = "Correct HDRs", alpha = "Correct HDRs")
# theme_bw()
plogrec_isomap_kde <-
  frank1 %>%
  ggplot(aes(x = fxy, y = fxy_hdr)) +
  # geom_rect(xmin = 0, xmax = rec[2], ymin = 0, ymax = rec[2], fill = greys[2], col = greys[2]) +
  # geom_rect(xmin = rec[2], xmax = rec[3], ymin = rec[2], ymax = rec[3], fill = greys[3], col = greys[3]) + # third blocks
  geom_rect(xmin = rec[1], xmax = rec[2], ymin = rec[1], ymax = rec[2], fill = greys[2], col = greys[2]) + # small blocks
  geom_rect(xmin = 0, xmax = rec[1], ymin = 0, ymax = rec[1], fill = greys[1], col = greys[1]) +
  # geom_point()+ # if not color HDRs
  geom_point(aes(col = col_hdr, shape = kde_false, size = kde_false, alpha = kde_false)) + # if color HDRs
  scale_color_manual(
    name = "HDRs",
    breaks = levels(frank1$col), # paste0(c(paste(sort(100 - prob)), ">99"), "%"),
    values = rev(c(RColorBrewer::brewer.pal(length(prob)+1, "YlOrRd")[-1], "#000000"))[1:length(levels)],
    labels = paste0(levels(frank1$col), "%")
  ) +
  scale_x_continuous(trans = "log10", limits = c(1, 2000), breaks = c(0, 20, 200, 1000, 2000), guide = guide_axis(n.dodge = 1)) + 
  scale_y_continuous(trans = "log10", limits = c(1, 2000), breaks = c(0, 20, 200, 1000, 2000)) + 
  coord_fixed() +
  scale_shape_manual(values = c(`TRUE` = 20, `FALSE` = 17 )) +
  scale_size_manual(values = c(`TRUE` = 1, `FALSE` = 2 )) +
  scale_alpha_manual(values = c(`TRUE` = .5, `FALSE` = 1) ) + 
  labs(x = "", y = "", shape = "Correct HDRs", size = "Correct HDRs", alpha = "Correct HDRs")

plotrec <- ((plogrec_isomap_vkde + labs(y = "Estimated density rank") + plogrec_isomap_kde) + 
              plot_layout(guides = 'collect')) %>%
  add_global_label(Xlab = "True density rank")
plotrec
ggsave(filename = paste0("paper/figures/", mapping, N,"_densityrank_comparison_", xname, "_radius", radius, "_r", format(r, decimal.mark = "_"), "_logrank_rec_colprob_smallblocks_crossfalse.png"), plotrec, width = 12, height = 5, dpi = 300)


# # Not run
# ##---------------------------------------------------------------------
# # color the points with distance to y=x (Improvement: rotate 45 degrees for better comparison)
# # log scale of x and y axis, with colors indicating the absolute rank error weighted by the sum of true and estimated ranks
# p_isomap_logrank_distcol <- frank %>% 
#   select(fxy, fxy_vkde, fxy_hdr) %>%
#   # mutate(outlier = case_when((fxy <= noutliers) ~ "True outliers",
#   #                            (fxy_dckde <= noutliers) ~ "DC-KDE outliers",
#   #                            (fxy_kde <= noutliers) ~ "KDE outliers",
#   #                            TRUE ~ "Not outliers")) %>% 
#   pivot_longer(cols = -c(fxy), names_to = "kde", values_to = "densities") %>%
#   mutate(fcol = abs((fxy - densities) / (fxy + densities) ) ) %>% # colors for weighted rank error
#   ggplot(aes(fxy, densities)) + # , col = densities
#   geom_point(aes(col = fcol)) + 
#   facet_grid(~ factor(kde, levels = c("fxy_vkde", "fxy_hdr")),
#              labeller = as_labeller(c(
#                `fxy_vkde` = paste("DC-KDE correlation =", format(cor(frank$fxy, frank$fxy_vkde), digits = 3)),
#                `fxy_hdr` = paste("KDE correlation =", format(cor(frank$fxy, frank$fxy_hdr), digits = 3))))
#   ) + 
#   scale_x_continuous(trans = "log10", limits = c(1, 2000), breaks = c(1, 10, 100, 1000, 2000)) +
#   scale_y_continuous(trans = "log10", limits = c(1, 2000), breaks = c(1, 10, 100, 1000, 2000)) +
#   labs(x = "True density rank", y = "Estimated density rank") + 
#   scale_color_viridis(option = "D", direction = -1) +
#   guides(color = "none")
# # ggsave(filename = paste0("paper/figures/", mapping, N,"_densityrank_comparison_isomap_radius", radius, "_r", format(r, decimal.mark = "_"), "_rankdistcol.png"), p_isomap_logrank_distcol, width = 12, height = 6, dpi = 300)
# 
# 
# # Not run
# ##---------------------------------------------------------------------
# # color the points with 4 kenrels
# pf_vkde <- frank %>% 
#   ggplot(aes(x = fxy, y = fxy_vkde, col = factor(preswissroll$label), shape = factor(preswissroll$label))) + 
#   geom_point() + 
#   labs(x = "", y = "", color = "Kernels", shape = "Kernels", title = paste("DC-KDE correlation =", format(cor(f$fxy_vkde, f$fxy, method = "spearman"), digits = 3)) ) +
#   scale_x_continuous(trans = "log10", limits = c(1, 2000), breaks = c(1, 10, 100, 1000, 2000)) +
#   scale_y_continuous(trans = "log10", limits = c(1, 2000), breaks = c(1, 10, 100, 1000, 2000))
# pf_hdr <- frank %>% 
#   ggplot(aes(x = fxy, y = fxy_hdr, col = factor(preswissroll$label), shape = factor(preswissroll$label))) + 
#   geom_point() + 
#   labs(x = "", y = "", color = "Kernels", shape = "Kernels", title = paste("KDE correlation =", format(cor(f$fxy_hdr, f$fxy, method = "spearman"), digits = 3))) +
#   scale_x_continuous(trans = "log10", limits = c(1, 2000), breaks = c(1, 10, 100, 1000, 2000)) +
#   scale_y_continuous(trans = "log10", limits = c(1, 2000), breaks = c(1, 10, 100, 1000, 2000))
# ((pf_vkde + labs(x = "True density rank", y = "Estimated density rank")) + (pf_hdr + labs(x = "True density rank"))) + 
#   coord_fixed() + 
#   plot_layout(guides = 'collect') & 
#   theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5, face = "bold"))
# # ggsave(paste0("paper/figures/", mapping, N,"_densityrank_comparison_isomap_radius", radius, "_r", format(r, decimal.mark = "_"), "_logrank_col4kernels.png"), width = 10, height = 6, dpi = 300)


# # DC-KDE (NOT RUN)
# frank1 <- frank %>% 
#   mutate(col_vkde = factor(100 - p_isomap$p[["plot_env"]][["region"]], 
#                            levels = c( "0", paste(sort(prob)))),
#          col_truehdr = factor(100 - p_meta$p[["plot_env"]][["region"]],
#                               levels = c( "0", paste(sort(prob))))
#          # col_truehdr = quantile(1:N, p_meta$p[["plot_env"]][["region"]]/100)
#   )
# rec <- quantile(1:N, c(sort(100-prob), 100)/100) %>% log10()
# greys <- scales::grey_pal()(10)[-(1:4)]
# # rec_data <- tibble(
# #   greys = scales::grey_pal()(10)[-(1:4)],
# #   xmin = 0,
# #   xmax = quantile(1:N, c(sort(100-prob), 100)/100)
# # )
# plogrec_isomap_dckde <- frank1 %>%
#   ggplot(aes(x = fxy, y = fxy_vkde)) + 
#   # geom_rect(fill=scales::grey_pal()(length(levels(frank1$col_truehdr))+5)[as.numeric(frank1$col_truehdr)+5],xmin = 0, xmax = quantile(1:N, c(sort(100-prob), 100)/100), ymin = -Inf,ymax = Inf, alpha = .1) +
#   # geom_rect(data = rec_data, aes(x = xmin, xmax = xmax, ymin = xmin, ymax = xmax, fill = greys, col = greys))
#   geom_rect(xmin = 0, xmax = rec[6], ymin = -Inf, ymax = rec[6], fill = greys[6], col = greys[6]) +
#   geom_rect(xmin = 0, xmax = rec[5], ymin = -Inf, ymax = rec[5], fill = greys[5], col = greys[5]) +
#   geom_rect(xmin = 0, xmax = rec[4], ymin = -Inf, ymax = rec[4], fill = greys[4], col = greys[4]) +
#   geom_rect(xmin = 0, xmax = rec[3], ymin = -Inf, ymax = rec[3], fill = greys[3], col = greys[3]) +
#   geom_rect(xmin = 0, xmax = rec[2], ymin = -Inf, ymax = rec[2], fill = greys[2], col = greys[2]) +
#   geom_rect(xmin = 0, xmax = rec[1], ymin = -Inf, ymax = rec[1], fill = greys[1], col = greys[1]) +
#   geom_point(aes(col = col_vkde)) +
#   scale_x_continuous(trans = "log10", limits = c(1, 2000), breaks = quantile(1:N, c(sort(100-prob)) / 100), guide = guide_axis(n.dodge = 2)) + 
#   scale_y_continuous(trans = "log10", limits = c(1, 2000), breaks = quantile(1:N, c(sort(100-prob)) / 100)) + # breaks = c(1, 10, 100, 1000, 2000)) +
#   # geom_point(aes(y = fxy_hdr, col = fxy_hdr)) +
#   # scale_color_viridis(option = "A", direction = -1) 
#   ggplot2::scale_color_manual(
#     name = "HDRs",
#     breaks = c(paste(sort(100 - prob)), "100"),
#     values = c(RColorBrewer::brewer.pal(length(prob)+1, "YlOrRd")[-1], "#000000")
#   ) + 
#   coord_fixed() +
#   labs(x = "True density quantiles", y = "Estimated density quantiles")
# # theme(axis.text.x = element_text(angle = 90, hjust = 1))






## -----------------------------------------------------------------------------
# 6.4 Plot of outlier rank to show that the ranks do not change across different ML methods
# change i to put each ML method on x-axis
i <- 1
xname <- c("isomap", "lle", "le", "tsne", "umap")[i]
mllabel <- c("ISOMAP", "LLE", "Laplacian Eigenmaps", "t-SNE", "UMAP")
p_data <- frank %>% 
  rename(fxy_vkde_isomap = fxy_vkde, 
         fxy_hdr_isomap = fxy_hdr) %>%
  # `colnames<-`(gsub("fxy_", "", colnames(.))) %>%
  select(-fxy) %>% 
  pivot_longer(-contains(paste0("_", xname)), names_sep = "_", names_to = c("f", "dckde",
                                                                            # pivot_longer(everything(), names_sep = "_", names_to = c("dckde",
                                                                            "kde")) %>% 
  pivot_longer(contains(xname), names_sep = "_", names_to = c("f_x",
                                                              "dckde_x", 
                                                              "kde_x"), 
               values_to = "x") %>% 
  filter(dckde == dckde_x)  %>% 
  pivot_wider(names_from = "kde_x", values_from = "x") %>% 
  select(kde, any_of(xname), value, dckde) %>% 
  mutate(dckde = factor(dckde, levels = c("vkde", "hdr"), labels = c("DC-KDE", "KDE")),
         kde = factor(kde, labels = mllabel[-i])) 
cor_data <- p_data %>%
  group_by(kde, dckde) %>% 
  summarise(cor = cor(value, !!sym(xname), method = "spearman"), 
            # x = quantile(!!sym(xname), 0.9), 
            # y = quantile(!!sym(xname), 0.2),
            .groups = "drop")

p_isovs4ml <- ggplot(p_data, aes(x = !!sym(xname), y = value)) +
  geom_point() +
  # geom_label(data = cor_data, aes(x=x, y = y, label = paste("cor: ", format(cor, digits = 2)))) +
  facet_grid(kde ~ dckde) + 
  labs(x = mllabel[i], y = "")
# ggsave(paste0("paper/figures/", mapping, N,"_density_compare_", xname, "vs4ml_radius", radius, "_r", format(r, decimal.mark = "_"), "_rank.png"), p_isovs4ml, width = 8, height = 10, dpi = 300)



# Not run
# frank
# prank_dckde_isolle <- frank %>% 
#   ggplot(aes(x = fxy_vkde, y = fxy_vkde_lle)) +  # , col = factor(preswissroll$label), shape = factor(preswissroll$label)
#   geom_point() + 
#   labs(x = "", y = "", subtitle = paste("Rank correlation =", format(cor(f$fxy_vkde, f$fxy_vkde_lle, method = "spearman"), digits = 3)) ) #  color = "Kernels", shape = "Kernels",
# prank_hdr_isolle <- frank %>% 
#   ggplot(aes(x = fxy_hdr, y = fxy_hdr_lle)) + 
#   geom_point() + 
#   labs(x = "", y = "", subtitle = paste("Rank correlation =", format(cor(f$fxy_hdr, f$fxy_hdr_lle, method = "spearman"), digits = 3))) #+
#   # scale_y_continuous(n.breaks = 5)
# 
# prank_dckde_isole <- frank %>% 
#   ggplot(aes(x = fxy_vkde, y = fxy_vkde_le)) + 
#   geom_point() + 
#   labs(x = "", y = "", subtitle = paste("Rank correlation =", format(cor(f$fxy_vkde, f$fxy_vkde_le, method = "spearman"), digits = 3)) )
# prank_hdr_isole <- frank %>% 
#   ggplot(aes(x = fxy_hdr, y = fxy_hdr_le)) + 
#   geom_point() + 
#   labs(x = "", y = "", subtitle = paste("Rank correlation =", format(cor(f$fxy_hdr, f$fxy_hdr_le, method = "spearman"), digits = 3))) # +
#   # scale_y_continuous(limits = c(0, max(f$fxy_hdr)), n.breaks = 5)
# 
# prank_dckde_isotsne <- frank %>% 
#   ggplot(aes(x = fxy_vkde, y = fxy_vkde_tsne)) + 
#   geom_point() + 
#   labs(x = "", y = "", subtitle = paste("Rank correlation =", format(cor(f$fxy_vkde, f$fxy_vkde_tsne, method = "spearman"), digits = 3)) )
# prank_hdr_isotsne <- frank %>% 
#   ggplot(aes(x = fxy_hdr, y = fxy_hdr_tsne)) + 
#   geom_point() + 
#   labs(x = "", y = "", subtitle = paste("Rank correlation =", format(cor(f$fxy_hdr, f$fxy_hdr_tsne, method = "spearman"), digits = 3, nsmall = 3))) # +
#   # scale_y_continuous(limits = c(0, max(f$fxy_hdr)), n.breaks = 5)
# 
# prank_dckde_isoumap <- frank %>% 
#   ggplot(aes(x = fxy_vkde, y = fxy_vkde_umap)) + 
#   geom_point() + 
#   labs(x = "", y = "", subtitle = paste("Rank correlation =", format(cor(f$fxy_vkde, f$fxy_vkde_umap, method = "spearman"), digits = 3)) )
# prank_hdr_isoumap <- frank %>% 
#   ggplot(aes(x = fxy_hdr, y = fxy_hdr_umap)) + 
#   geom_point() + 
#   labs(x = "", y = "", subtitle = paste("Rank correlation =", format(cor(f$fxy_hdr, f$fxy_hdr_umap, method = "spearman"), digits = 3))) # +
#   # scale_y_continuous(limits = c(0, max(f$fxy_hdr)), n.breaks = 5)
# 
# (prank_dckde_isolle + prank_hdr_isolle)
# 
# prank_isomap <- (
#     (prank_dckde_isolle + labs(y = "LLE", title = "DC-KDE") + (prank_hdr_isolle + labs(title = "KDE") + theme(axis.text.y = element_blank()))) / 
#     (prank_dckde_isole + ylab("Laplacian Eigenmaps") + prank_hdr_isole + theme(axis.text.y = element_blank())) / 
#     (prank_dckde_isotsne + ylab("t-SNE") + prank_hdr_isotsne + theme(axis.text.y = element_blank())) / 
#     (prank_dckde_isoumap + ylab("UMAP") + xlab("ISOMAP") + (prank_hdr_isoumap + xlab("ISOMAP") + theme(axis.text.y = element_blank()))) 
# ) +
#   plot_layout(guides = 'collect') & 
#   theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.4, face = "bold"))
# prank_isomap
# # ggsave(paste0("paper/figures/", mapping, N,"_density_compare_isomapvs4ml_radius", radius, "_r", format(r, decimal.mark = "_"), "_rank.png"), prank_isomap, width = 8, height = 10, dpi = 300)


# Not run
# Plot all combinations of dc-kde / kde in one plot
# library(GGally)
# ggpairs(frank[,(1:5)*2],
#         columnLabels = c("ISOMAP", "LLE", "Laplacian Eigenmaps", "t-SNE", "UMAP"),
#         upper = list(continuous = wrap('cor', size = 3)),
#         diag = "blank",
#         title = "Rank comparison of DC-KDE for different manifold learning methods")
# ggpairs(frank[,(1:5)*2+1],
#         columnLabels = c("ISOMAP", "LLE", "Laplacian Eigenmaps", "t-SNE", "UMAP"),
#         upper = list(continuous = wrap('cor', size = 3)),
#         diag = "blank",
#         title = "Rank comparison of KDE for different manifold learning methods")


# TODO:
# Spec and sensitivity for correctly categorized outliers, over a small area in the top left (Done using log10 instead)
# Change plot of true density colors to 5 levels (Done)
# Table of rank correlation (Done)
# Table of embedding quality
