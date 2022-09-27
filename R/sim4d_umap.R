## -----------------------------------------------------------------------------
## This script is for manifold learning and kernel density estimation (KDE and DC-KDE)
## -----------------------------------------------------------------------------
rm(list = ls())
library(tidyverse)
library(dimRed)
library(hdrcde)
library(ggforce)
library(ks)
library(patchwork)
Jmisc::sourceAll(here::here("R/sources"))
set.seed(1234)
r <- 1
# # r <- as.numeric(commandArgs()[[6]])
load("data/simdata_100d_4dmanifold_N10000_trueden_k200.rda")

paste("Start at:", Sys.time())
# ----parameters----------------------------------------------------------------
x <- train
y <- NULL # without pre-computed embedding
N <- nrow(x)
s <- 4 # embedded in 4-D
k <- N / 50
method <- "annUMAP"
annmethod <- "kdtree"
distance <- "euclidean"
treetype <- "kd"
searchtype <- "radius" # change searchtype for radius search based on `radius`, or KNN search based on `k`
radius <- 10 # the bandwidth parameter, \sqrt(\elsilon), as in algorithm. Note that the radius need to be changed for different datasets, not to increase k

gridsize <- 20
opt.method <- c("AMISE", "MEAN", "SCALED")[3] # 3 ONLY FOR NOW, no scaling for Rn
riem.scale <- 1 # tune parameter

## ----umap-------------------------------------------------------------------
metric_umap <- metricML(x, fn = y, s = s, k = k, radius = radius, method = method, invert.h = TRUE, eps = 0,
                          annmethod = annmethod, distance = distance, treetype = treetype,
                          searchtype = searchtype
)
fn <- metric_umap$embedding
Rn <- metric_umap$rmetric
adj_matrix <- metric_umap$adj_matrix
E1 <- fn[,1]; E2 <- fn[,2]

print("True density summary statistics")
summary(trueden)

## ----Fixed bandwidth KDE with ks::kde-----------------------------------------
tictoc::tic()
fixden_umap <- vkde(x = fn, h = NULL, vh = NULL, r = r, gridsize = gridsize, eval.points = fn, opt.method = opt.method, riem.scale = riem.scale, adj_matrix = adj_matrix)
tictoc::toc()
(H <- fixden_umap$H)
print("KDE summary statistics")
summary(fixden_umap$estimate)


## ----DCKDE--------------------------------------------------------------------
# r <- sqrt(mean(apply(Rn, 3, det)))
# r <- 1 # bandwidth parameter in vkde() # r = 0.1, too small, many true densities are estimated as 0; r = 2, too large, many zero true densities are estimated as non-zeros; r = 0.5 ~ 0.8, 
tictoc::tic()
fumap <- vkde(x = fn, h = NULL, vh = Rn, r = r, gridsize = gridsize, eval.points = fn, opt.method = opt.method, riem.scale = riem.scale, adj_matrix = adj_matrix)
tictoc::toc()
print("DC-KDE summary statistics")
summary(fumap$estimate)


# ----Compare with true density-------------------------------------------------
cormethod <- "spearman"
(cors <- c(
  cor(trueden, fumap$estimate, method = cormethod),
  cor(trueden, fixden_umap$estimate, method = cormethod)
))
mean(sqrt((trueden - fumap$estimate)^2)) # smaller MSE
mean(sqrt((trueden - fixden_umap$estimate)^2))

noutliers <- 100
sum(head(order(fumap$estimate), noutliers) %in% head(order(trueden), noutliers))
sum(head(order(fixden_umap$estimate), noutliers) %in% head(order(trueden), noutliers))

f <- tibble(fxy = trueden, # true densities
            fxy_dckde = fumap$estimate, fxy_hdr = fixden_umap$estimate
) %>% summarise_all(rank) # Comment this part if plotting the density estimates instead of ranks
f
# Plot ranks instead of densities
p <- f %>% 
  pivot_longer(cols = -1, names_to = "kde", values_to = "densities") %>% 
  ggplot(aes(x = fxy, y = densities)) + 
  geom_point() + 
  facet_grid(~ kde, labelle = as_labeller(c("fxy_dckde" = paste("DC-KDE correlation", round(cor(f$fxy_dckde, f$fxy, method = "spearman"), 3)), 
                                            "fxy_hdr" = paste("KDE correlation", round(cor(f$fxy_hdr, f$fxy, method = "spearman"), 3))))
  ) + 
  labs(x = "True density rank", y = "Estimated density rank") +
  scale_y_continuous(n.breaks = 6)
# p
ggsave(paste0("figures/compareden_5d_N", N, "_", method, "_r", format(r, decimal.mark = "_"), ".png"), p, width = 10, height = 6, dpi = 300)


save(method, fixden_umap, fumap, train, trueden, cors,
     file = paste0("data/compareden_4d_N", N, "_", method, "_radius", radius, "_r", format(r, decimal.mark = "_"), "_annumap.rda"))
save(metric_umap, file = paste0("data/metric_umap_4d_N10000_radius", radius, ".rda"))

paste("End at:", Sys.time())
