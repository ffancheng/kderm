###--------------------------------------------------------
## 1. Gather and plot the results from 4 scenarios in HPC results
###--------------------------------------------------------
rm(list = ls())
library(tidyverse)
library(viridis)
N <- 10000
radius <- 8
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

# getOption("digits")
# options(digits = 7)
f <- tibble(fxy = fxy,
            fxy_vkde_isomap = fisomap$estimate, fxy_hdr_isomap = fixden_isomap$estimate,
            fxy_vkde_lle = flle$estimate, fxy_hdr_lle = fixden_lle$estimate,
            fxy_vkde_le = fle$estimate, fxy_hdr_le = fixden_le$estimate,
            fxy_vkde_umap = fumap$estimate, fxy_hdr_umap = fixden_umap$estimate
)
# f[which(rank(f$fxy_vkde_isomap) == 33.5),]$fxy_hdr_isomap # There are ties in dc-kde results, leading to rank=33.5 for the top 67 anomalies
# f[which(rank(f$fxy_vkde_isomap) == 33.5),]$fxy # Their true densities are not exactly the same, but are all below 1e-07, so basically 0

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
## 2.2 Table for correctly identified outliers percentage
###--------------------------------------------------------
prob <- c(1, 10, 50, 95, 99)
hdrbreaks <- c(0, quantile(1:N, c(sort(100 - prob) / 100)), N)
hdrlabels <- c(">99", paste0(100 - prob)) #, "%", not working for latex, \\%
frankhdr <- f %>% 
  summarise_all(.fun = list(~ rank(x = ., ties.method = "average"))) %>% # using average ties
  # mutate(outliers = factor(label, levels = c("TRUE", "FALSE")),
  #        dist2center = dist2center) %>%
  `colnames<-`(gsub("fxy_", "", colnames(.))) %>% 
  pivot_longer(cols = -c(fxy), names_to = "kde", values_to = "densities") %>% 
  separate(kde, c("KDE", "ML"), sep = "_") %>% 
  mutate(KDE = factor(KDE, levels = c("vkde", "hdr"), labels = c("DC-KDE", "KDE")),
         ML = factor(ML, levels = methods, labels = c("ISOMAP", "LLE", "Laplacian Eigenmaps", "UMAP"))
  ) %>% 
  mutate(truehdr = cut(fxy, 
                       breaks = hdrbreaks, 
                       labels = hdrlabels),
         esthdr = cut(densities, 
                      breaks = hdrbreaks, 
                      labels = hdrlabels)
         ) %>% 
  mutate(correcthdr = (truehdr == esthdr))
hdrtable <- frankhdr %>% 
  group_by(KDE, ML, truehdr) %>% 
  summarise(numhdr = sum(correcthdr) / n(), .groups = "drop") %>% 
  pivot_wider(names_from = c(KDE, ML), values_from = numhdr) %>% 
  column_to_rownames(var = "truehdr") %>%
  `colnames<-` (NULL)
hdrmax <- (as.matrix(hdrtable[, 1:4]) - as.matrix(hdrtable[, 5:8])) >= 0
hdrmax <- cbind(hdrmax, !hdrmax)[, c(1, 5, 2, 6, 3, 7, 4, 8)] %>% `rownames<-` (NULL)
hdrtable <- hdrtable[, c(1, 5, 2, 6, 3, 7, 4, 8)]
hdrtable %>% 
  head(n = 2) %>% 
  `rownames<-` (c(">99\\% HDR", "99\\% HDR")) %>% 
  kableExtra::kbl(caption = "Percentage comparison of correct highest density regions in density estimation of four manifold learning embeddings.", booktabs = TRUE, digits = 3, escape = FALSE, align = "c",
                  # col.names = c(rep(c("DC-KDE", "KDE"), 4))
                  ) %>%
  kable_classic("striped") %>% 
  kable_styling(latex_options = "scale_down") %>%
  kable_paper(full_width = TRUE) %>%
  add_header_above(c("HDRs" = 1, rep(c("DC-KDE" = 1, "KDE" = 1), 4))) %>%
  add_header_above(c(" " = 1, "ISOMAP" = 2, "LLE" = 2, "Laplacian Eigenmaps" = 2, "UMAP" = 2)) %>% 
  # column_spec(1, bold = TRUE) %>% 
  # add_indent(1:6) %>% 
  column_spec(2, bold = hdrmax[,1]) %>% 
  column_spec(3, bold = hdrmax[,2]) %>% 
  column_spec(4, bold = hdrmax[,3]) %>% 
  column_spec(5, bold = hdrmax[,4]) %>% 
  column_spec(6, bold = hdrmax[,5]) %>% 
  column_spec(7, bold = hdrmax[,6]) %>% 
  column_spec(8, bold = hdrmax[,7]) %>% 
  column_spec(2*(1:4), border_left = T, background = grey(0.9))

# save(frankhdr, hdrtable, hdrmax, file = paste0("paper/figures/hdrtable_", "4d_N", N, "_4ml_radius", 10, "_r", 1, ".rda"))


###--------------------------------------------------------
## 3. Rank comparison plot, highlight outliers
###--------------------------------------------------------
frank <- f %>% 
  summarise_all(.fun = list(~ rank(x = ., ties.method = "average"))) %>% # random tie method for better visualization (not horizontal lines in >99% region)
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
  scale_size_manual(values = c(`TRUE` = 2, `FALSE` = 1), labels = c(2,1)) +
  scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 20), labels = c(2,1)) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = .2), labels = c(2,1) ) +
  scale_color_viridis(option = "inferno", direction = -1) +
  labs(x = "True density rank", y = "Estimated density rank", col = "Distance", size = "Mixture components", shape = "Mixture components", alpha = "Mixture components") +
  scale_x_log10() + 
  scale_y_log10() + 
  theme(legend.position = "bottom")
p_rank
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
#   labs(x = "", y = "Laplacian Eigenmaps", color = "Distance", shape = "Kernels", title = "") +
#   scale_y_continuous(limits = c(0, max(f$fxy_vkde)), n.breaks = 6)
# pf_hdr <- f %>%
#   ggplot(aes(x = fxy, y = fxy_hdr, col = dist2center, shape = label)) +
#   # scale_color_viridis(option = "A", direction = -1) +
#   geom_point() +
#   labs(x = "", y = "", color = "Distance", shape = "Kernels", title = "") +
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
