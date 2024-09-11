###--------------------------------------------------------
## 1. Gather and plot the results from 4 scenarios in HPC results: 4-D manifold embedding in 100-D space simulation
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
# label[which(rank(f$fxy_vkde_isomap) == 33.5)] # lowest density points are from both mixture components

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
# ----KDE with ks::kde on 100-d data instead of 4-d embedding-------------------
library(weird)
library(tictoc)

# From Rob
n <- 100
d <- 3
y <- matrix(rnorm(n*d), ncol = d)
H <- kde_bandwidth(y)
fhat <- kde(y, H)


# Handcode KDE
library(mvtnorm)
multivariate_kde <- function(x, eval_points, H = NULL) {
  n <- nrow(x)
  d <- ncol(x)
  m <- nrow(eval_points)
  
  # If H is not provided, use Scott's rule
  if (is.null(H)) {
    H <- diag(d) * n^(-1/(d+4))
  }
  
  # Ensure H is symmetric and positive definite
  H <- (H + t(H)) / 2
  eigen_values <- eigen(H)$values
  if (any(eigen_values <= 0)) {
    stop("Bandwidth matrix H must be positive definite")
  }
  
  # Compute the inverse and determinant of H
  H_inv <- solve(H)
  H_det <- det(H)
  
  # Compute the density estimate
  density <- numeric(m)
  for (i in 1:m) {
    for (j in 1:n) {
      diff <- eval_points[i,] - x[j,]
      density[i] <- density[i] + exp(-0.5 * t(diff) %*% H_inv %*% diff)
    }
  }
  
  # Normalize the density estimate
  density <- density / (n * (2*pi)^(d/2) * sqrt(H_det))
  
  return(density)
}
# # Generate some sample data
# set.seed(123)
# n <- 1000
# d <- 3
# x <- matrix(rnorm(n*d), ncol=d)
# # Define points to evaluate the density
# eval_points <- matrix(seq(-3, 3, length.out=20), ncol=d, nrow=20, byrow=TRUE)
# # Estimate the density
# density_estimates <- multivariate_kde(x, eval_points)
# print(density_estimates)



# Simplify the function with H=diag(d)
multivariate_kde_I <- function(x, eval_points) {
  n <- nrow(x)
  d <- ncol(x)
  m <- nrow(eval_points)
  
  # Use identity matrix for H
  H <- diag(d)
  
  # Compute the inverse and determinant of H (which are simple for identity matrix)
  H_inv <- H  # Inverse of identity matrix is itself
  H_det <- 1  # Determinant of identity matrix is 1
  
  # Compute the density estimate
  density <- numeric(m)
  for (i in 1:m) {
    for (j in 1:n) {
      diff <- eval_points[i,] - x[j,]
      density[i] <- density[i] + exp(-0.5 * sum(diff^2))  # Simplified due to identity H
    }
  }
  
  # Normalize the density estimate
  density <- density / (n * (2*pi)^(d/2))  # Simplified due to H_det = 1
  
  return(density)
}
all.equal(
  multivariate_kde_I(x, eval_points),
  multivariate_kde(x, eval_points, H = diag(d))
)



# set H=diag(d)
# X: true manifold, train: 100d data for dimRed, N, p, den: density of 5d semishepre, trueden: true density of 4d manifold

# KDE on true manifold X
dim(X) # 10000*4
dim(train) # 10000*100
tic()
# H <- kde_bandwidth(X)
fhat <- multivariate_kde_I(x = train, eval_points = train)
toc() # 15min
summary(fhat)

# set H as Identity matrix for 100d data
H1 <- diag(ncol(train)) #+ diag(1e-8, ncol(train))
tic()
fhat1 <- kde(x = train, H = H1, binned = FALSE, eval.points = train)
toc() # 20min
summary(fhat1$estimate)

all.equal(fhat, fhat1$estimate) # TRUE

# True manifold density
summary(trueden)
# fxy <- trueden
# dencor <- function(x) cor(x$estimate, fxy, method = "spearman")
dencor(fhat1) # 0.9884194
cor(fhat, trueden, method = "s") # 0.9884194


# Scale H
optimal_h <- function(X, method = "scott") {
  n <- nrow(X)  # number of observations
  d <- ncol(X)  # number of dimensions
  
  if (method == "scott") {
    # Scott's rule
    h <- n^(-1 / (d + 4))
  } else if (method == "silverman") {
    # Silverman's rule (modified for multivariate case)
    h <- (4 / (d + 2))^(1 / (d + 4)) * n^(-1 / (d + 4))
  } else if (method == "normal_reference") {
    # Normal reference rule
    sigma <- apply(X, 2, sd)
    h <- (4 / (d + 2))^(1 / (d + 4)) * n^(-1 / (d + 4)) * mean(sigma)
  } else {
    stop("Unknown method. Please choose 'scott', 'silverman', or 'normal_reference'.")
  }
  
  return(h)
}
X = train
h_scott <- optimal_h(X, method = "scott")
h_silverman <- optimal_h(X, method = "silverman")
h_normal <- optimal_h(X, method = "normal_reference")


# Now use H = h^2 * I, use ks::kde()
H <- h_silverman ^ 2 * diag(ncol(train))
tic()
fhat0 <- kde(x = train, H = H, binned = FALSE, eval.points = train)
toc() # 20min
summary(fhat0$estimate)
cor(fhat0$estimate, trueden, method = "s") # 0.9884194



summary(trueden)
dencor(fhat0) # 0.9882927



# ----Compare with true density and kde, dckde----------------------------------
# Table 2
col.max <- apply(cors, 2, max)
max.match <- rep(col.max, each = 2) == cors
data.frame(cors) %>%
  mutate(across(where(is.numeric), ~format(round(., 3), nsmall = 3))) %>%
  kableExtra::kbl(caption = "Correlation between true density and estimated density for four manifold learning embeddings of the 100-D data. DC-KDE has a slightly lower rank correlation for ISOMAP and LLE above 0.96, but significantly outperfoms KDE for Laplacian Eigenmaps and UMAP when more distortion is induced in dimension reduction.", booktabs = TRUE, digits = 3, escape = FALSE) %>%
  # kable_styling(latex_options = "scale_down") %>% # automatically adjust the table to page width
  kable_paper(full_width = FALSE) %>%
  # column_spec(1, width = "6cm") %>%
  column_spec(2, bold = max.match[,1]) %>%
  column_spec(3, bold = max.match[,2]) %>%
  column_spec(4, bold = max.match[,3]) %>%
  column_spec(5, bold = max.match[,4])


# Calculate the dencor value
kde_no_embedding_value <- dencor(fhat1)
cors %>%
  rbind(`KDE no embedding` = rep(kde_no_embedding_value, ncol(.))) %>%
  as.data.frame() %>%
  mutate(across(everything(), ~format(round(as.numeric(.), 3), nsmall = 3))) %>%
  kbl(caption = "Correlation between true density and estimated density for four manifold learning embeddings of the 100-D data. DC-KDE has a slightly lower rank correlation for ISOMAP and LLE above 0.96, but significantly outperfoms KDE for Laplacian Eigenmaps and UMAP when more distortion is induced in dimension reduction. KDE no embedding is added for comparison.", 
      booktabs = TRUE, 
      digits = 3, 
      escape = FALSE) %>%
  kable_paper(full_width = FALSE) %>%
  column_spec(1, bold = c(max.match[,1], FALSE)) %>%
  column_spec(2, bold = c(max.match[,2], FALSE)) %>%
  column_spec(3, bold = c(max.match[,3], FALSE)) %>%
  column_spec(4, bold = c(max.match[,4], FALSE)) %>%
  row_spec(nrow(updated_cors), hline_after = TRUE)




# Table 3
frankhdr1 <- f %>%
  mutate(fhat = fhat1$estimate) |> 
  summarise_all(.fun = list(~ rank(x = ., ties.method = "average"))) %>% # using average ties
  # mutate(outliers = factor(label, levels = c("TRUE", "FALSE")),
  #        dist2center = dist2center) %>%
  `colnames<-`(gsub("fxy_", "", colnames(.))) %>% 
  pivot_longer(cols = -c(fxy), names_to = "kde", values_to = "densities") %>% 
  separate(kde, c("KDE", "ML"), sep = "_") %>% 
  mutate(KDE = factor(KDE, levels = c("vkde", "hdr", "fhat"), labels = c("DC-KDE", "KDE", "KDE no embedding")),
         ML = factor(ML, levels = c(methods, ""), labels = c("ISOMAP", "LLE", "Laplacian Eigenmaps", "UMAP", ""))
  ) %>% 
  mutate(truehdr = cut(fxy, 
                       breaks = hdrbreaks, 
                       labels = hdrlabels),
         esthdr = cut(densities, 
                      breaks = hdrbreaks, 
                      labels = hdrlabels)
  ) %>% 
  mutate(correcthdr = (truehdr == esthdr))
hdrtable <- 
  frankhdr1 %>% 
  group_by(KDE, ML, truehdr) %>% 
  summarise(numhdr = sum(correcthdr) / n(), .groups = "drop") %>% 
  pivot_wider(names_from = c(KDE, ML), values_from = numhdr) %>% 
  column_to_rownames(var = "truehdr") %>%
  `colnames<-` (NULL)
hdrmax <- (as.matrix(hdrtable[, 1:4]) - as.matrix(hdrtable[, 5:8])) >= 0
hdrmax <- cbind(hdrmax, !hdrmax)[, c(1, 5, 2, 6, 3, 7, 4, 8)] %>% `rownames<-` (NULL)
hdrtable <- hdrtable[, c(1, 5, 2, 6, 3, 7, 4, 8, 9)]
hdrtable %>% 
  head(n = 2) %>% 
  `rownames<-` (c(">99\\% HDR", "99\\% HDR")) %>% 
  kableExtra::kbl(caption = "Percentage comparison of correct highest density regions in density estimation of four manifold learning embeddings.", booktabs = TRUE, digits = 3, escape = FALSE, align = "c",
                  # col.names = c(rep(c("DC-KDE", "KDE"), 4))
  ) %>%
  kable_classic("striped") %>% 
  kable_styling(latex_options = "scale_down") %>%
  kable_paper(full_width = TRUE) %>%
  add_header_above(c("HDRs" = 1, rep(c("DC-KDE" = 1, "KDE" = 1), 4), "KDE" = 1)) %>%
  add_header_above(c(" " = 1, "ISOMAP" = 2, "LLE" = 2, "Laplacian Eigenmaps" = 2, "UMAP" = 2, "No embedding" = 1)) %>% 
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




# Figure 9, row 1 only
frank <- f %>% 
  select(fxy, ends_with("isomap")) %>%
  mutate(fhat = fhat1$estimate) %>%
  summarise_all(.fun = list(~ rank(x = ., ties.method = "average"))) %>% # random tie method for better visualization (not horizontal lines in >99% region)
  mutate(outliers = factor(label, levels = c("TRUE", "FALSE")),
         dist2center = dist2center)

frank_long <- frank %>%
  `colnames<-`(gsub("fxy_", "", colnames(.))) %>% 
  pivot_longer(cols = -c(fxy, outliers, dist2center), names_to = "kde", values_to = "densities") %>% 
  separate(kde, c("KDE", "ML"), sep = "_") %>% 
  mutate(KDE = factor(KDE, levels = c("vkde", "hdr", "fhat"), labels = c("DC-KDE", "KDE", "KDE no embedding")),
         # ML = factor(ML, levels = methods, labels = c("ISOMAP", "LLE", "Laplacian Eigenmaps", "UMAP"))
  ) |> 
  select(-ML) |> 
  distinct()


p_rank <- 
  frank_long %>%
  ggplot(aes(fxy, densities, col = dist2center)) +
  geom_rect(xmin = 0, xmax = log10(N / 100), ymin = 0, ymax = log10(N / 100), fill = "#B7B7B7", col = "#B7B7B7") +
  geom_point(aes(shape = outliers, size = outliers, alpha = outliers)) +
  facet_grid( ~ KDE) +
  scale_size_manual(values = c(`TRUE` = 2, `FALSE` = 1), labels = c(2,1)) +
  scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 20), labels = c(2,1)) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = .2), labels = c(2,1) ) +
  scale_color_viridis(option = "inferno", direction = -1) +
  labs(x = "True density rank", y = "Estimated density rank", col = "Distance", size = "Mixture components", shape = "Mixture components", alpha = "Mixture components") +
  scale_x_log10() + 
  scale_y_log10() + 
  theme(legend.position = "bottom")
p_rank





# Plot histogram of density estimates

frank <- f %>% 
  select(fxy, ends_with("isomap")) %>%
  mutate(fhat = fhat1$estimate) %>%
  # summarise_all(.fun = list(~ rank(x = ., ties.method = "average"))) %>% # random tie method for better visualization (not horizontal lines in >99% region)
  mutate(outliers = factor(label, levels = c("TRUE", "FALSE")),
         dist2center = dist2center)

frank_long <- frank %>%
  `colnames<-`(gsub("fxy_", "", colnames(.))) %>% 
  pivot_longer(cols = -c(outliers, dist2center), names_to = "kde", values_to = "densities") %>% 
  separate(kde, c("KDE", "ML"), sep = "_") %>% 
  mutate(KDE = factor(KDE, levels = c("fxy", "vkde", "hdr", "fhat"), labels = c("True density", "DC-KDE", "KDE on embedding", "KDE on 100-D data")),
         # ML = factor(ML, levels = methods, labels = c("ISOMAP", "LLE", "Laplacian Eigenmaps", "UMAP"))
  ) |> 
  select(-ML) |> 
  distinct()

# frank_long |> 
#   ggplot(aes(x = KDE, y = densities, color = KDE)) +
#   geom_point() +
#   geom_violin() +
#   # ylim(0, 0.1) +
#   facet_wrap(~ KDE, scales = "free_y", ncol = 1) +
#   theme(legend.position = "none")


# Put the violin plot in the review paper
library(viridis)
# library(hrbrthemes)  # For a sleek theme

frank_long |> 
  ggplot(aes(x = KDE, y = densities, fill = KDE, color = KDE)) +
  geom_violin(alpha = 0.7, scale = "width") +
  # geom_boxplot(width = 0.1, alpha = 0.5, fatten = NULL) +
  # geom_jitter(alpha = 0.1, width = 0.1, height = 0) +
  scale_y_continuous(n.breaks = 5) +
  facet_wrap(~ KDE, scales = "free_y", nrow = 4) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  scale_color_viridis(discrete = TRUE, option = "D") +
  labs(title = "Comparison of density estimates for the semi-hypersphere example",
       subtitle = "Distribution of densities for different KDE methods",
       x = "KDE Methods",
       y = "Densities") +
  theme_ipsum() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_blank(), # element_text(face = "bold", size = 12, hjust = 0),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray50"),
    axis.text.x = element_text(face = "bold"),
    panel.spacing = unit(1, "lines")
  )
ggsave("jcgs/revision/kde_compare_violin.jpeg", dpi = 300, width = 8, height = 10)





## -----------------------------------------------------------------------------











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
