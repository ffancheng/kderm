###--------------------------------------------------------
## 1. Gather and plot the results from 4 scenarios in HPC results
###--------------------------------------------------------
rm(list = ls())
# first load all .rda files
library(tidyverse)
library(viridis)
N <- 10000
radius <- 10
r <- 1
# method <- c("isomap", "lle", "le", "tsne", "umap")[1]
load("data/simdata_100d_4dmanifold_N10000_trueden_k100.rda")
fxy <- trueden
label <- as.factor(c(rep(FALSE, 0.99*N), rep(TRUE, 0.01*N))) # where the points are manual outliers
dist2center <- sqrt(radius^2 - semisphere[,5]^2) # Euclidean distance to the center O
methods <- c("isomap", "lle", "le", "umap")
# load HPC results for four methods
files <- list.files(path = "data/truedenk100", pattern = "*.rda", full.names = TRUE)
lapply(files, load, .GlobalEnv)


###--------------------------------------------------------
## 2. Table for density correlations
###--------------------------------------------------------
fxy <- trueden
dencor <- function(x) cor(x$estimate, fxy, method = "spearman")
# dencor(p_lle)
# dencor(p_hdr_lle)
cors <- cbind(
  c(dencor(fisomap), dencor(fixden_isomap)),
  c(dencor(flle), dencor(fixden_lle)),
  c(dencor(fle), dencor(fixden_le)),
  c(dencor(fumap), dencor(fixden_umap))
) 
rownames(cors) <- c("DC-KDE", "KDE")
colnames(cors) <- c("ISOMAP", "LLE", "Laplacian Eigenmaps", "UMAP")
cors %>% 
  kableExtra::kbl(caption = "Correlation between true density and estimated density for four manifold learning embeddings.", booktabs = TRUE, digits = 3, escape = FALSE) %>%
  kable_styling(latex_options = "scale_down") %>%
  kable_paper(full_width = TRUE) 
# save(cors, file = paste0("paper/figures/CorrelationTable_", "4d_N", N, "_4ml_radius", 10, "_r", 1, ".rda"))


###--------------------------------------------------------
## 2. Table for correctlly identified outliers percentage
###--------------------------------------------------------




###--------------------------------------------------------
## 3. Rank comparison plot, highlight outliers
###--------------------------------------------------------
f <- tibble(fxy = fxy,
            fxy_vkde_isomap = fisomap$estimate, fxy_hdr_isomap = fixden_isomap$estimate,
            fxy_vkde_lle = flle$estimate, fxy_hdr_lle = fixden_lle$estimate,
            fxy_vkde_le = fle$estimate, fxy_hdr_le = fixden_le$estimate,
            fxy_vkde_umap = fumap$estimate, fxy_hdr_umap = fixden_umap$estimate
            )
frank <- f %>% 
  summarise_all(.fun = list(~ rank(x = ., ties.method = "random"))) %>% 
  mutate(outliers = factor(label, levels = c("TRUE", "FALSE")),
         dist2center = dist2center)

frank_long <- frank %>%
  `colnames<-`(gsub("fxy_", "", colnames(.))) %>% 
  pivot_longer(cols = -c(fxy, outliers, dist2center), names_to = "kde", values_to = "densities") %>% 
  separate(kde, c("KDE", "ML"), sep = "_") %>% 
  mutate(KDE = factor(KDE, levels = c("vkde", "hdr"), labels = c("DC-KDE", "KDE")),
         ML = factor(ML, levels = methods, labels = c("ISOMAP", "LLE", "Laplacian Eigenmaps", "UMAP"))
  )
p_rank <- frank_long %>%
  ggplot(aes(fxy, densities, col = dist2center)) +
  geom_rect(xmin = 0, xmax = log10(N / 100), ymin = 0, ymax = log10(N / 100), fill = "#B7B7B7", col = "#B7B7B7") +
  geom_point(aes(shape = outliers, size = outliers, alpha = outliers)) +
  facet_grid(ML ~ KDE) +
  scale_size_manual(values = c(`TRUE` = 2, `FALSE` = 1)) +
  scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 20)) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = .2) ) + 
  # scale_color_manual(
  #   name = "HDRs",
  #   breaks = levels(frank$col),
  #   values = rev(c(RColorBrewer::brewer.pal(2, "YlOrRd")[-1], "#000000"))[1:2],
  #   labels = paste0(levels(frank$col), "%")
  # )
  scale_color_viridis(option = "inferno", direction = -1) +
  labs(x = "True density rank", y = "Estimated density rank", col = "Distance", size = "Outliers", shape = "Outliers", alpha = "Outliers") +
  scale_x_log10() + 
  scale_y_log10() + 
  theme(legend.position = "bottom")
# ggsave(paste0("paper/figures/", "sim4d", N, "_density_comparison_4ml_radius", 10, "_k200_rankdensity_circleoutlier_with1rec.png"), p_rank, width = 8, height = 10, dpi = 300)

# Not run
# ## scatterplot to compare f_xy for ISOMAP
# f <- tibble(fxy = fxy, fxy_vkde = fisomap$estimate, fxy_hdr = fixden_isomap$estimate) %>% summarise_all(rank)
# f1 <- tibble(fxy = fxy, fxy_vkde = flle$estimate, fxy_hdr = fixden_lle$estimate) %>% summarise_all(rank)
# f2 <- tibble(fxy = fxy, fxy_vkde = fle$estimate, fxy_hdr = fixden_le$estimate) %>% summarise_all(rank)
# f3 <- tibble(fxy = fxy, fxy_vkde = fumap$estimate, fxy_hdr = fixden_umap$estimate) %>% summarise_all(rank)
# 
# f=f
# pf_vkde <- f %>%
#   ggplot(aes(x = fxy, y = fxy_vkde, col = dist2center, shape = label)) +
#   # scale_color_viridis(option = "A", direction = -1) +
#   geom_point() +
#   labs(x = "", y = "ISOMAP", color = "Distance", shape = "Kernels", title = "DC-KDE") +
#   scale_y_continuous(limits = c(0, max(f$fxy_vkde)), n.breaks = 6)
# pf_hdr <- f %>%
#   ggplot(aes(x = fxy, y = fxy_hdr, col = dist2center, shape = label)) +
#   geom_point() +
#   labs(x = "", y = "", color = "Distance", shape = "Kernels", title = "KDE") +
#   scale_y_continuous(limits = c(0, max(f$fxy_hdr)), n.breaks = 6)
# p <- pf_vkde + pf_hdr
# 
# f=f1
# pf_vkde <- f %>%
#   ggplot(aes(x = fxy, y = fxy_vkde, col = dist2center, shape = label)) +
#   # scale_color_viridis(option = "A", direction = -1) +
#   geom_point() +
#   labs(x = "", y = "LLE", color = "Distance", shape = "Kernels", title = "") +
#   scale_y_continuous(limits = c(0, max(f$fxy_vkde)), n.breaks = 6)
# pf_hdr <- f %>%
#   ggplot(aes(x = fxy, y = fxy_hdr, col = dist2center, shape = label)) +
#   # scale_color_viridis(option = "A", direction = -1) +
#   geom_point() +
#   labs(x = "", y = "", color = "Distance", shape = "Kernels", title = "") +
#   scale_y_continuous(limits = c(0, max(f$fxy_hdr)), n.breaks = 6)
# p1 <- pf_vkde + pf_hdr
# 
# f=f2
# pf_vkde <- f %>%
#   ggplot(aes(x = fxy, y = fxy_vkde, col = dist2center, shape = label)) +
#   # scale_color_viridis(option = "A", direction = -1) +
#   geom_point() +
#   labs(x = "True density", y = "Laplacian Eigenmaps", color = "Distance", shape = "Kernels", title = "") +
#   scale_y_continuous(limits = c(0, max(f$fxy_vkde)), n.breaks = 6)
# pf_hdr <- f %>%
#   ggplot(aes(x = fxy, y = fxy_hdr, col = dist2center, shape = label)) +
#   # scale_color_viridis(option = "A", direction = -1) +
#   geom_point() +
#   labs(x = "True density", y = "", color = "Distance", shape = "Kernels", title = "") +
#   scale_y_continuous(limits = c(0, max(f$fxy_hdr)), n.breaks = 6)
# p2 <- pf_vkde + pf_hdr
# 
# f=f3
# pf_vkde <- f %>%
#   ggplot(aes(x = fxy, y = fxy_vkde, col = dist2center, shape = label)) +
#   geom_point() +
#   labs(x = "", y = "UMAP", color = "Distance", shape = "Kernels", title = "") +
#   scale_y_continuous(limits = c(0, max(f$fxy_vkde)), n.breaks = 6)
# pf_hdr <- f %>%
#   ggplot(aes(x = fxy, y = fxy_hdr, col = dist2center, shape = label)) +
#   geom_point() +
#   labs(x = "", y = "", color = "Distance", shape = "Kernels", title = "") +
#   scale_y_continuous(limits = c(0, max(f$fxy_hdr)), n.breaks = 6)
# p3 <- pf_vkde + pf_hdr
# 
# result <- (p/p1/p2/p3) + plot_layout(guides = 'collect') &
#   theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5, face = "bold")) &
#   scale_color_viridis(option = "A", direction = -1)
# gt <- patchwork::patchworkGrob(result)
# gt <- gridExtra::grid.arrange(gt, left = "Estimated density", bottom = "True density")
# gt
