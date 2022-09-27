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
method <- "annIsomap"
annmethod <- "kdtree"
distance <- "euclidean"
treetype <- "kd"
searchtype <- "radius" # change searchtype for radius search based on `radius`, or KNN search based on `k`
radius <- 10 # the bandwidth parameter, \sqrt(\elsilon), as in algorithm. Note that the radius need to be changed for different datasets, not to increase k

gridsize <- 20
opt.method <- c("AMISE", "MEAN", "SCALED")[3] # 3 ONLY FOR NOW, no scaling for Rn
riem.scale <- 1 # tune parameter

## ----isomap-------------------------------------------------------------------
metric_isomap <- metricML(x, fn = y, s = s, k = k, radius = radius, method = method, invert.h = TRUE, eps = 0,
                          annmethod = annmethod, distance = distance, treetype = treetype,
                          searchtype = searchtype
)
fn <- metric_isomap$embedding
Rn <- metric_isomap$rmetric
adj_matrix <- metric_isomap$adj_matrix
E1 <- fn[,1]; E2 <- fn[,2]

print("True density summary statistics")
summary(trueden)

## ----Fixed bandwidth KDE with ks::kde-----------------------------------------
tictoc::tic()
fixden_isomap <- vkde(x = fn, h = NULL, vh = NULL, r = r, gridsize = gridsize, eval.points = fn, opt.method = opt.method, riem.scale = riem.scale, adj_matrix = adj_matrix)
tictoc::toc()
(H <- fixden_isomap$H)
print("KDE summary statistics")
summary(fixden_isomap$estimate)


## ----DCKDE--------------------------------------------------------------------
# r <- sqrt(mean(apply(Rn, 3, det)))
# r <- 1 # bandwidth parameter in vkde() # r = 0.1, too small, many true densities are estimated as 0; r = 2, too large, many zero true densities are estimated as non-zeros; r = 0.5 ~ 0.8, 
tictoc::tic()
fisomap <- vkde(x = fn, h = NULL, vh = Rn, r = r, gridsize = gridsize, eval.points = fn, opt.method = opt.method, riem.scale = riem.scale, adj_matrix = adj_matrix)
tictoc::toc()
print("DC-KDE summary statistics")
summary(fisomap$estimate)


# ----Compare with true density-------------------------------------------------
cormethod <- "spearman"
(cors <- c(
  cor(trueden, fisomap$estimate, method = cormethod),
  cor(trueden, fixden_isomap$estimate, method = cormethod)
))
mean(sqrt((trueden - fisomap$estimate)^2)) # smaller MSE
mean(sqrt((trueden - fixden_isomap$estimate)^2))

noutliers <- 100
sum(head(order(fisomap$estimate), noutliers) %in% head(order(trueden), noutliers))
sum(head(order(fixden_isomap$estimate), noutliers) %in% head(order(trueden), noutliers))

f <- tibble(fxy = trueden, # true densities
            fxy_dckde = fisomap$estimate, fxy_hdr = fixden_isomap$estimate
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


save(method, fixden_isomap, fisomap, train, trueden, cors,
     file = paste0("data/compareden_4d_N", N, "_", method, "_radius", radius, "_r", format(r, decimal.mark = "_"), "_annIsomap.rda"))
save(metric_isomap, file = paste0("data/metric_isomap_4d_N10000_radius", radius, ".rda"))

paste("End at:", Sys.time())
