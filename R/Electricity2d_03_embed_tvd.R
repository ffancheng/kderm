rm(list = ls())
library(tidyverse)
library(dimRed)
library(dtplyr)
library(data.table)
library(reticulate)
library(viridis)
# library(ks)
library(hdrcde)
library(igraph)
library(matrixcalc)
library(ggforce)
library(patchwork)
Jmisc::sourceAll("R/sources")
# scen <- 1
scen <- as.numeric(commandArgs()[[6]])
r <- 180 #1 
# quantile(nn2res$nn.dists)
# 0%      25%      50%      75%     100% 
# 0.0000 179.6016 187.3877 194.8748 277.0549 

load("data/half_count_ratio_3639id336tow.rda")
ids <- spdemand$id
train <- spdemand[, !"id"]
rm(spdemand)

paste("Start at:", Sys.time())
# ----parameters----------------------------------------------------------------
# x <- train
y <- NULL # without pre-computed embedding
N <- nrow(train)
s <- 2 # embedded in 2-D for visualization purpose, 6-D for estimated intrinsic dimension
k <- 100 # N / 50
method <- c("annIsomap", "annLLE", "annLaplacianEigenmaps", "anntSNE", "annUMAP")[scen]
annmethod <- "kdtree"
distance <- c("euclidean", "manhattan")[2] # "manhattan" for all households
treetype <- "kd"
searchtype <- "priority" # "radius" # change searchtype for radius search based on `radius`, or KNN search based on `k`
radius <- 20 # the bandwidth parameter, \sqrt(\elsilon), as in algorithm. Note that the radius need to be changed for different datasets, not to increase k

gridsize <- 20
opt.method <- c("AMISE", "MEAN", "SCALED")[3] # 3 ONLY FOR NOW, no scaling for Rn
riem.scale <- 1 # tune parameter
prob <- c(1, 10, 50, 95, 99)
noutliers <- 20

## ----embedding----------------------------------------------------------------
paste("Embedding began at:", Sys.time())
pars <- list(knn = k,
             radius = radius,
             eps = 0, # for true NN 
             ndim = s,
             get_geod = FALSE,
             annmethod = "kdtree", 
             nt = 50, 
             search.k = 50,
             nlinks = 16, 
             ef.construction = 500, 
             distance = distance,
             treetype = treetype,
             searchtype = searchtype)
isomapnn <- embed(train, .method = method, knn = pars$k, radius = pars$radius,
      annmethod = pars$annmethod,
      eps = pars$eps,
      nt = pars$nt,
      nlinks = pars$nlinks, ef.construction = pars$ef.construction,
      distance = pars$distance,
      treetype = pars$treetype,
      searchtype = pars$searchtype,
      ndim = pars$ndim,
      .mute = c("output"))
# ann_table_isomapnn <- calc_ann_table(e = isomapnn, nn.idx = truenn)
# dr_quality(X=e@org.data, Y=e@data@data, K=e@pars$knn, nn.idx = nn.idx)$quality
fn <- isomapnn@data@data

# # Use existing results for embedding
# load(paste0("data/metric_", method, "_electricity_2d_radius20.rda"))
# fn <- metric_isomap$embedding


## ---learnmetric---------------------------------------------------------------
paste("Learn metric began at:", Sys.time())
# if(file.exists("data/metric_isomap_4d_N10000_radius10.rda")){
#   load("data/metric_isomap_4d_N10000_radius10.rda")
# } else {
  metric_isomap <- metricML(x = train, fn = fn, bandwidth = 4, #!!!
                            s = s, k = k, radius = radius, method = method, invert.h = TRUE, eps = 0,
                            annmethod = annmethod, distance = distance, treetype = treetype,
                            searchtype = searchtype
  )
# }
# fn <- metric_isomap$embedding
save(metric_isomap, file = paste0("data/metric_", method, "_electricity_2d_radius", radius, "_k", k, "_searchtype", searchtype,  "_r", format(r, decimal.mark = "_"), ".rda"))
Rn <- metric_isomap$rmetric
adj_matrix <- metric_isomap$adj_matrix
E1 <- fn[,1]; E2 <- fn[,2]
head(Rn)
paste("nn2res dimension", str(metric_isomap$nn2res))
head(metric_isomap$nn2res$nn.idx)

## ----Fixed bandwidth KDE with ks::kde-----------------------------------------
paste("KDE began at:", Sys.time())
tictoc::tic()
fixden_isomap <- vkde(x = fn, h = NULL, vh = NULL, r = r, gridsize = gridsize, eval.points = fn, opt.method = opt.method, riem.scale = riem.scale, adj_matrix = adj_matrix)
tictoc::toc()
print("KDE summary statistics with vkde NULL h")
summary(fixden_isomap$estimate)
p_hdr_isomap <- plot_outlier(x = list(embedding = fn, Rn = NULL), gridsize = gridsize, prob = prob, riem.scale = riem.scale, f = fixden_isomap, ell.size = 0, label = ids)
# Or plot using hdrcde
# tictoc::tic()
# p_hdr_isomap <- hdrscatterplot_new(fn[,1], fn[,2], kde.package = "ks", levels = prob, noutliers = noutliers, label = ids)
# tictoc::toc()
# print("KDE summary statistics with hdrcde")
# summary(p_hdr_isomap$densities)

## ----DCKDE--------------------------------------------------------------------
paste("DC-KDE began at:", Sys.time())
tictoc::tic()
fisomap <- vkde(x = fn, h = NULL, vh = Rn, r = r, gridsize = gridsize, eval.points = fn, opt.method = opt.method, riem.scale = riem.scale, adj_matrix = adj_matrix) ### TODO: NAN?????
tictoc::toc()
print("DC-KDE summary statistics")
summary(fisomap$estimate)
p_isomap <- plot_outlier(x = metric_isomap, gridsize = gridsize, prob = prob, riem.scale = riem.scale, f = fisomap, ell.size = 0, label = ids)

# p_hdr_isomap$p + p_isomap$p

save(method, fixden_isomap, fisomap, p_hdr_isomap, p_isomap,
     file = paste0("data/compareden_electricity_2d_N", N, "_", method, "_radius", radius, "_k", k, "_searchtype", searchtype, "_r", format(r, decimal.mark = "_"), ".rda"))
# save(metric_isomap, file = paste0("data/metric_", method, "_electricity_2d_radius", radius, "_k", k, "_searchtype", searchtype,  "_r", format(r, decimal.mark = "_"), ".rda"))

paste("End at:", Sys.time())

# Not run (no true density)
# ----Compare with true density-------------------------------------------------
# cormethod <- "spearman"
# (cors <- c(
#   cor(trueden, fisomap$estimate, method = cormethod),
#   cor(trueden, fixden_isomap$estimate, method = cormethod)
# ))
# mean(sqrt((trueden - fisomap$estimate)^2)) # smaller MSE
# mean(sqrt((trueden - fixden_isomap$estimate)^2))
# 
# noutliers <- 100
# sum(head(order(fisomap$estimate), noutliers) %in% head(order(trueden), noutliers))
# sum(head(order(fixden_isomap$estimate), noutliers) %in% head(order(trueden), noutliers))
# 
# f <- tibble(fxy = trueden, # true densities
#             fxy_dckde = fisomap$estimate, fxy_hdr = fixden_isomap$estimate
# ) %>% summarise_all(rank) # Comment this part if plotting the density estimates instead of ranks
# f
# # Plot ranks instead of densities
# p <- f %>% 
#   pivot_longer(cols = -1, names_to = "kde", values_to = "densities") %>% 
#   ggplot(aes(x = fxy, y = densities)) + 
#   geom_point() + 
#   facet_grid(~ kde, labelle = as_labeller(c("fxy_dckde" = paste("DC-KDE correlation", round(cor(f$fxy_dckde, f$fxy, method = "spearman"), 3)), 
#                                             "fxy_hdr" = paste("KDE correlation", round(cor(f$fxy_hdr, f$fxy, method = "spearman"), 3))))
#   ) + 
#   labs(x = "True density rank", y = "Estimated density rank") +
#   scale_y_continuous(n.breaks = 6)
# # p
# ggsave(paste0("figures/compareden_electricity_2d_N", N, "_", method, "_r", format(r, decimal.mark = "_"), ".png"), p, width = 10, height = 6, dpi = 300)
