###--------------------------------------------------------
## Gather and plot the results from 4 senarios in HPC results
###--------------------------------------------------------
# first load all .rda files
library(viridis)
N <- 10000
radius <- 10
r <- 1
# method <- c("isomap", "lle", "le", "tsne", "umap")[1]
load("~/git/kderm/data/simdata_3d_N2000_01_trueden_k100.rda")
fxy <- trueden
label <- as.factor(c(rep(1, 0.99*N), rep(2, 0.01*N)))
dist2center <- sqrt(radius^2 - semisphere[,5]^2)
methods <- c("isomap", "lle", "le", "umap")
## scatterplot to compare f_xy for ISOMAP
f <- tibble(fxy = fxy, fxy_vkde = fisomap$estimate, fxy_hdr = fixden_isomap$estimate) #%>% rank()
f1 <- tibble(fxy = fxy, fxy_vkde = flle$estimate, fxy_hdr = fixden_lle$estimate) #%>% rank()
f2 <- tibble(fxy = fxy, fxy_vkde = fle$estimate, fxy_hdr = fixden_le$estimate) #%>% rank()
f3 <- tibble(fxy = fxy, fxy_vkde = fumap$estimate, fxy_hdr = fixden_umap$estimate) #%>% rank()

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
ggsave(paste0("paper/figures/", "sim4d", N, "_density_comparison_4ml_riem", format(r, decimal.mark = "_"), "k200_rankdensity.png"), gt, width = 8, height = 10, dpi = 300)



# Table for density correlations
fxy <- trueden
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
