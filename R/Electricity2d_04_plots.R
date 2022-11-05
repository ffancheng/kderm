# This script contains the code for plotting for smart meter data
rm(list=ls())
library(data.table)
library(dtplyr)
library(tidyverse)
library(dimRed)
library(ggplot2)
library(patchwork)
library(hdrcde)
library(igraph)
library(MASS)
library(ggraph)
library(pryr)
library(colorspace)
library(caret)
library(tictoc)
library(ggforce)
library(kableExtra)
Jmisc::sourceAll(here::here("R/sources")) 
set.seed(1)

# scen <- 1
# scen <- as.numeric(commandArgs()[[6]])

###------------------------------------------------------
nid <- 3639
ntow <- 336
len <- 100
filename <- paste0(nid, "id_", ntow, "tow_", len, "length")
load(paste0('data/half_count_ratio_3639id336tow.rda')) # 3639*(336*201+1)
train <- spdemand[, !"id"]
ids <- spdemand$id
N <- nrow(train)
folder <- "data/"
distance <- c("euclidean", "manhattan")[2]
###------------------------------------------------------
# method <- c("annIsomap", "annLLE", "annLaplacianEigenmaps", "anntSNE", "annUMAP")[scen]
N <- nrow(train)
s <- 2 # embedded in 2-D for visualization purpose
k <- 100 # N / 50
radius <- 20 # the bandwidth parameter, \sqrt(\elsilon), as in algorithm. Note that the radius need to be changed for different datasets, not to increase k
searchtype <- "priority" # "radius" # change searchtype for radius search based on `radius`, or KNN search based on `k`
annmethod <- "kdtree"
distance <- c("euclidean", "manhattan")[2] # "manhattan" for all households
treetype <- "kd"

r <- 180
gridsize <- 20
opt.method <- c("AMISE", "MEAN", "SCALED")[3] # 3 ONLY FOR NOW, no scaling for Rn
riem.scale <- 1 
prob <- c(1, 10, 50, 95, 99)
noutliers <- 20

par <- 0 # eps = 0
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

# load("data/truenn_3639id_tvd.rda")
# load("data/truenn_tvd_3639id_336tow_100length_radius20_k100_searchtypepriority.rda")



###-------------------------------------------------
### Load 5 Ml results
###-------------------------------------------------
ml <- c("annIsomap", "annLLE", "annLaplacianEigenmaps", "anntSNE", "annUMAP")

method <- ml[2]
load( paste0("data/compareden_electricity_2d_N", N, "_", method, "_radius", radius, "_k", k, "_searchtype", searchtype, "_r", format(r, decimal.mark = "_"), ".rda"))
flle <- fisomap
fixden_lle <- fixden_isomap
p_hdr_lle <- p_hdr_isomap
p_lle <- p_isomap
load(paste0("data/metric_", method, "_electricity_2d_radius", radius, "_k", k, "_searchtype", searchtype,  "_r", format(r, decimal.mark = "_"), ".rda"))
metric_lle <- metric_isomap

method <- ml[3]
load( paste0("data/compareden_electricity_2d_N", N, "_", method, "_radius", radius, "_k", k, "_searchtype", searchtype, "_r", format(r, decimal.mark = "_"), ".rda"))
fle <- fisomap
fixden_le <- fixden_isomap
p_hdr_le <- p_hdr_isomap
p_le <- p_isomap
load(paste0("data/metric_", method, "_electricity_2d_radius", radius, "_k", k, "_searchtype", searchtype,  "_r", format(r, decimal.mark = "_"), ".rda"))
metric_le <- metric_isomap

method <- ml[4]
load( paste0("data/compareden_electricity_2d_N", N, "_", method, "_radius", radius, "_k", k, "_searchtype", searchtype, "_r", format(r, decimal.mark = "_"), ".rda"))
ftsne <- fisomap
fixden_tsne <- fixden_isomap
p_hdr_tsne <- p_hdr_isomap
p_tsne <- p_isomap
load(paste0("data/metric_", method, "_electricity_2d_radius", radius, "_k", k, "_searchtype", searchtype,  "_r", format(r, decimal.mark = "_"), ".rda"))
metric_tsne <- metric_isomap

method <- ml[5]
load( paste0("data/compareden_electricity_2d_N", N, "_", method, "_radius", radius, "_k", k, "_searchtype", searchtype, "_r", format(r, decimal.mark = "_"), ".rda"))
fumap <- fisomap
fixden_umap <- fixden_isomap
p_hdr_umap <- p_hdr_isomap
p_umap <- p_isomap
load(paste0("data/metric_", method, "_electricity_2d_radius", radius, "_k", k, "_searchtype", searchtype,  "_r", format(r, decimal.mark = "_"), ".rda"))
metric_umap <- metric_isomap

method <- ml[1]
load(paste0("data/compareden_electricity_2d_N", N, "_", method, "_radius", radius, "_k", k, "_searchtype", searchtype, "_r", format(r, decimal.mark = "_"), ".rda"))
# fisomap <- fisomap
# fixden_isomap <- fixden_isomap
# p_hdr_isomap <- p_hdr_isomap
# p_isomap <- p_isomap
load(paste0("data/metric_", method, "_electricity_2d_radius", radius, "_k", k, "_searchtype", searchtype,  "_r", format(r, decimal.mark = "_"), ".rda"))
# metric_isomap <- metric_isomap

cat("Load .rda finished at: "); Sys.time()



###-------------------------------------------------
### Highest density regions plots
###-------------------------------------------------
# 2D embedding plot with highest density regions
nolabs <- labs(x = "", y = "")
noytext <-  theme(axis.text.y = element_blank()) 
phdr5ml <- (
  ((
       (p_isomap$p + labs(x = "", y = "ISOMAP", title = "DC-KDE")) |
       (p_hdr_isomap$p + labs(x = "", y = "", title = "KDE") + noytext ) 
  ) & scale_y_continuous(breaks = c(-250, 0, 250))) /
    (( 
         (p_lle$p + labs(x = "", y = "LLE")) |
         (p_hdr_lle$p + nolabs + noytext) 
    ) & scale_y_continuous(breaks = c(0, 5, 10)) & scale_x_continuous(n.breaks = 4)) /
    (( 
         (p_le$p + labs(x = "", y = "Laplacian Eigenmaps")) |
         (p_hdr_le$p + nolabs + noytext) 
    ) & scale_x_continuous(breaks = c(-.025, 0, .025)) & scale_y_continuous(breaks = c(-.05, 0, .05), limits = c(-.05, .05)) ) /
    (( 
         (p_tsne$p + labs(x = "", y = "t-SNE")) |
         (p_hdr_tsne$p + nolabs + noytext) 
    ) & scale_y_continuous(breaks = c(-20, 0, 20)) & scale_x_continuous(n.breaks = 5)) / 
    (( 
         (p_umap$p + labs(x = "", y = "UMAP")) |
         (p_hdr_umap$p + nolabs + noytext) 
    ) & scale_y_continuous(n.breaks = 4))
) + 
  plot_layout(guides = 'collect') &
  # guides(HDRs = guide_legend(nrow = 1)) &
  theme(legend.direction = "horizontal", legend.position = "bottom", 
        legend.box = "horizontal", axis.title.y = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"))
phdr5ml
ggsave(paste0("paper/figures/Electricity_2d_5levels_outliers_comparison_5ml_r", format(r, decimal.mark = "_"), ".png"), phdr5ml, width = 9, height = 12, dpi = 300)


###-------------------------------------------------
## Table for top 20 outlier IDs
###-------------------------------------------------
# ids
# p_isomap$outlier
# p_hdr_isomap$outlier
# p_hdr_isomap$densities # unordered
# Outliers
f <- tibble(
  fxy_vkde_isomap = ids[p_isomap$outlier], fxy_hdr_isomap = ids[p_hdr_isomap$outlier],
  fxy_vkde_lle = ids[p_lle$outlier], fxy_hdr_lle = ids[p_hdr_lle$outlier],
  fxy_vkde_le = ids[p_le$outlier], fxy_hdr_le = ids[p_hdr_le$outlier],
  fxy_vkde_tsne = ids[p_tsne$outlier], fxy_hdr_tsne = ids[p_hdr_tsne$outlier],
  fxy_vkde_umap = ids[p_umap$outlier], fxy_hdr_umap = ids[p_hdr_umap$outlier]
)
f
f20 <- f[,c(1,3,5,7,9, 2,4,6,8,10)] %>% 
  head(n = 20) %>% 
  rowid_to_column()
f20rep <- f20[,-1] %>% unlist() %>% table() %>% as.data.frame() %>% filter(Freq > 3) %>% arrange(desc(Freq))
freqmax <- apply(f20, 2, function(x) x %in% f20rep$.) # highlight all ids that appears >3 times
freqmax1 <- apply(f20, 2, function(x) x %in% f20rep$.[f20rep$Freq == 5])
freqmax2 <- apply(f20, 2, function(x) x %in% f20rep$.[f20rep$Freq == 4])
f20cols <- freqmax1*2 + freqmax2
f20cols <- f20cols %>% 
  as_tibble() %>%
  mutate_all(as.character) %>% 
  mutate_all(funs(str_replace(., "1", blues9[5]))) %>% 
  mutate_all(funs(str_replace(., "2", blues9[9]))) %>% 
  as.matrix()
tabletop20 <- f20 %>% 
  kableExtra::kbl(caption = "Top 20 anomalous household IDs using density estimates of different 2-dimensional embedding. The left five columns are anomalies by DC-KDE, while the right five columns are anomalies using KDE. The blue colored IDs are IDs that are detected over three times, with the darker one indicating five times and the lighter one indicating four times.", booktabs = TRUE, digits = 3, escape = FALSE, col.names = NULL, align=rep('c', 10)) %>%
  kable_styling(latex_options = "scale_down") %>%
  kable_paper(full_width = TRUE) %>%
  add_header_above(c(" " = 1, rep(c("ISOMAP" = 1, "LLE" = 1, "Laplacian Eigenmaps" = 1, "t-SNE" = 1, "UMAP" = 1), 2))) %>% 
  add_header_above(c(" " = 1, "DC-KDE" = 5, "KDE" = 5)) %>%
  # row_spec(1, bold = TRUE) %>% 
  column_spec(1, bold = FALSE) %>% 
  column_spec(2, bold = freqmax[,2], color = f20cols[,2]) %>% 
  column_spec(3, bold = freqmax[,3], color = f20cols[,3]) %>% 
  column_spec(4, bold = freqmax[,4], color = f20cols[,4]) %>% 
  column_spec(5, bold = freqmax[,5], color = f20cols[,5]) %>% 
  column_spec(6, bold = freqmax[,6], color = f20cols[,6]) %>% 
  column_spec(7, bold = freqmax[,7], color = f20cols[,7], border_left = TRUE) %>% 
  column_spec(8, bold = freqmax[,8], color = f20cols[,8]) %>% # background = grey(0.9)
  column_spec(9, bold = freqmax[,9], color = f20cols[,9]) %>% 
  column_spec(10, bold = freqmax[,10], color = f20cols[,10]) %>% 
  column_spec(11, bold = freqmax[,11], color = f20cols[,11])
tabletop20
save(f20, f20cols, freqmax, tabletop20, file = paste0("paper/figures/Electricity_2d_table_outlier_idcompare_5ml_radius", radius, "_r", format(r, decimal.mark = "_"), ".rda"))




## -----------------------------------------------------------------------------
# Plot of outlier rank to show that the ranks do not change across different ML methods
## -----------------------------------------------------------------------------
# Densities
fden <- tibble(
  fxy_vkde_isomap = p_isomap$densities, fxy_hdr_isomap = p_hdr_isomap$densities,
  fxy_vkde_lle = p_lle$densities, fxy_hdr_lle = p_hdr_lle$densities,
  fxy_vkde_le = p_le$densities, fxy_hdr_le = p_hdr_le$densities,
  fxy_vkde_tsne = p_tsne$densities, fxy_hdr_tsne = p_hdr_tsne$densities,
  fxy_vkde_umap = p_umap$densities, fxy_hdr_umap = p_hdr_umap$densities
)
fden
frank <- fden %>% summarise_all(rank)
# frank_long <- frank %>%
#   # mutate(outlier = case_when((fxy <= noutliers) ~ "True outliers",
#   #                            (fxy_dckde <= noutliers) ~ "DC-KDE outliers",
#   #                            (fxy_kde <= noutliers) ~ "KDE outliers",
#   #                            TRUE ~ "Not outliers")) %>%
#   `colnames<-`(gsub("fxy_", "", colnames(.))) %>% 
#   pivot_longer(everything(), names_to = "kde", values_to = "densities") %>% 
#   separate(kde, c("KDE", "ML"), sep = "_") %>% 
#   mutate(KDE = factor(KDE, levels = c("vkde", "hdr"), labels = c("DC-KDE", "KDE")),
#          ML = factor(ML, levels = c("isomap", "lle", "le", "tsne", "umap"), labels = c("ISOMAP", "LLE", "Laplacian Eigenmaps", "t-SNE", "UMAP"))
#   )
# gt <- frank_long %>%
#   ggplot(aes(fxy, densities)) +
#   geom_point() +
#   facet_grid(ML ~ KDE) +
#   # scale_x_log10() +
#   # scale_y_log10() +
#   labs(x = "True density rank", y = "Estimated density rank")
# gt


# change i to put each ML method on x-axis
# Use UMAP
i <- 5
xname <- c("isomap", "lle", "le", "tsne", "umap")[i]
mllabel <- c("ISOMAP", "LLE", "Laplacian Eigenmaps", "t-SNE", "UMAP")
p_data <- frank %>% 
  pivot_longer(-contains(paste0("_", xname)), names_sep = "_", names_to = c("f", "dckde",
                                                                            # pivot_longer(everything(), names_sep = "_", names_to = c("dckde",
                                                                            "kde")) %>% 
  pivot_longer(contains(xname), names_sep = "_", names_to = c("f_x",
                                                              "dckde_x", 
                                                              "kde_x"), 
               values_to = "x") %>% 
  filter(dckde == dckde_x)  %>% 
  pivot_wider(names_from = "kde_x", values_from = "x") %>% 
  dplyr::select(kde, any_of(xname), value, dckde) %>% 
  mutate(dckde = factor(dckde, levels = c("vkde", "hdr"), labels = c("DC-KDE", "KDE")),
         kde = factor(kde, labels = mllabel[-i])) 
# cor_data <- p_data %>%
#   group_by(kde, dckde) %>% 
#   summarise(cor = cor(value, !!sym(xname), method = "spearman"), 
#             # x = quantile(!!sym(xname), 0.9), 
#             # y = quantile(!!sym(xname), 0.2),
#             .groups = "drop")
p_isovs4ml <- ggplot(p_data, aes(x = !!sym(xname), y = value)) +
  geom_rect(xmin = 0, xmax = log10(20), ymin = 0, ymax = log10(20), fill = "#B7B7B7", col = "#B7B7B7") +
  geom_point() +
  # geom_label(data = cor_data, aes(x=x, y = y, label = paste("cor: ", format(cor, digits = 2)))) +
  facet_grid(kde ~ dckde) + 
  labs(x = mllabel[i], y = "") +
  scale_x_log10() + 
  scale_y_log10()
p_isovs4ml
ggsave(paste0("paper/figures/Electricity_2d_rankplot_", N,"id_", xname, "vs4ml_radius", radius, "_r", format(r, decimal.mark = "_"), ".png"), p_isovs4ml, width = 8, height = 10, dpi = 300)



###--------------------------------------------------------
## Table for correctly identified outliers percentage >99% & 99%
###--------------------------------------------------------
prob <- c(1, 10, 50, 95, 99)
hdrbreaks <- c(0, quantile(1:N, c(sort(100 - prob) / 100)), N)
hdrlabels <- c(">99", paste0(100 - prob)) #, "%", not working for latex, \\%
frankhdr <- p_data %>% 
  `colnames<-` (c("ML", "fxy", "densities", "KDE")) %>% 
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
table_outlierpercent <- 
hdrtable %>% 
  head(n = 2) %>% 
  `rownames<-` (c(">99\\% HDR", "99\\% HDR")) %>% 
  kableExtra::kbl(caption = "Percentage of the same outliers as UMAP in two highest density regions when estimating densities for four manifold learning embeddings.", booktabs = TRUE, digits = 3, escape = FALSE, align = "c",
                  # col.names = c(rep(c("DC-KDE", "KDE"), 4))
  ) %>%
  kable_classic("striped") %>% 
  kable_styling(latex_options = "scale_down") %>%
  kable_paper(full_width = TRUE) %>%
  add_header_above(c("HDRs" = 1, rep(c("DC-KDE" = 1, "KDE" = 1), 4))) %>%
  add_header_above(c(" " = 1, "ISOMAP" = 2, "LLE" = 2, "Laplacian Eigenmaps" = 2, "t-SNE" = 2)) %>% 
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
table_outlierpercent
save(frankhdr, hdrtable, hdrmax, table_outlierpercent, file = paste0("paper/figures/Electricity_2d_hdrtable_umapvs4ml_radius", radius, "_r", format(r, decimal.mark = "_"), ".rda"))



###--------------------------------------------------------
## Plot of three anomalous and 1 typical households
###--------------------------------------------------------
f
f20 <- f[,c(1,3,5,7,9, 2,4,6,8,10)] %>% 
  head(n = 20) %>% 
  rowid_to_column()
f20rep <- f20[,-1] %>% unlist() %>% table() %>% as.data.frame() %>% filter(Freq > 3) %>% arrange(desc(Freq))
f20rep$.[f20rep$Freq == 5]

ftail10 <- f[,c(1,3,5,7,9, 2,4,6,8,10)] %>% 
  tail(n = 20)
ftail10rep <- ftail10 %>% unlist() %>% table() %>% as.data.frame() %>% filter(Freq >= 2) %>% arrange(desc(Freq))
ftail10rep$.[ftail10rep$Freq == 4] # 2396 3683 4857 7043
# Plot 2396 as typical household, found 4 times using DC-KDE except for tsne

# Typical: 2396
ftail10rep$.[ftail10rep$Freq == 4]
# Anomalous: 3161 4546 7049
f20rep$.[f20rep$Freq == 5]
ids4plot <- c(ftail10rep$.[ftail10rep$Freq == 4], f20rep$.[f20rep$Freq == 5])

# # Run once
# # Plot the electricity demand distribution for four households
# load("data/DT.rda")
# head(DT)
# DT4id <- DT[id %in% ids4plot,]
# save(DT4id, file = "data/Electricity_2d_plot_DTcompare4id_1typical_3anomalous.rda")

###-------------------------------------------------
### Plot typical households
###-------------------------------------------------
load("data/Electricity_2d_plot_DTcompare4id_1typical_3anomalous.rda")
DT4id
summary(DT4id)
p4id <- DT4id %>%
  lazy_dt() %>% 
  mutate(id = factor(id, levels = ids4plot)) %>% 
  as.data.table() %>% 
  ggplot(aes(x = day, y = demand, group = id)) +
  geom_line() +
  facet_grid(id ~ .) +
  labs(x = "Days", y = "Demand (kWh)")
p4id
ggsave(paste0("paper/figures/Electricity_2d_compare4id_1typical_3anomalous.png"), p4id, width = 8, height = 8, dpi = 500)
save(p4id, file = "data/Electricityplot_compare4id.rda")

# # Not run
# # Plot the electricity demand distribution
# DT4id
# p4ids <- DT4id %>% 
#   as_tibble() %>% 
#   mutate(id = as.factor(id),
#          typical = ifelse(id == unique(DT4id$id)[1], TRUE, FALSE)) %>% 
#   ggplot(aes(x=demand, group=id, fill=typical)) + 
#   geom_histogram(
#     aes(y=..density..), # Histogram with density instead of count on y-axis
#     binwidth=.1,
#     colour="grey40", alpha=.2, # position="density"
#   ) + 
#   geom_density(alpha=.3) +  # Overlay with transparent density plot
#   scale_fill_manual(values = c("black", "orange")) +
#   facet_wrap(~id, ncol = 1, scales = "free_y", strip.position="right") + # ,labeller = as_labeller(c("4669"="ID 4669", "3243"="ID 3243", "1321"="ID 1321"))
#   labs(x = "Demand (kWh)", y = "Density") + 
#   theme(legend.position = "None")
# p4ids
# # wrap_plots(p, p4ids)
# pcomp <- cowplot::plot_grid(
#   p, p4ids,
#   align = "h", axis = "tb",
#   nrow = 1, rel_widths = c(2, 1)
# )
# pcomp
# ggsave(paste0("paper/figures/Electricity_2d_compare4id_histogram.png"), pcomp, width = 10, height = 6)
# # plotly::ggplotly(p, width = 600, height = 400)


# Not run
# ## KDE and DC-KDE raw time series plot
# load("data/electricityplot_compare3id.rda")
# p_dc <- p
# load("data/electricityplot_DTcompare3id.rda")
# ((p_dc + ggtitle("DC-KDE")) +
#     (p + ggtitle("KDE"))) + 
#   plot_layout(guides = "collect") &
#   theme(#legend.position = 'bottom', 
#     plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) 
# ggsave("paper/figures/dckde_2densities_raw_electricity.png")



###-------------------------------------------------
### Plot 4 households distributions
###-------------------------------------------------
# remotes::install_github("ffancheng/gghdr") # add nxmargin argument for number of boxplots
library(gghdr)
cde4id <- 
  DT4id %>% 
  as_tibble() %>% 
  mutate(tow = as.integer(tow),
         dow = str_replace_all(dow, c("1"="Monday", "2"="Tuesday", "3" = "Wednesday", "4" = "Thursday", "5" = "Friday", "6" = "Saturday", "7" = "Sunday")), 
         dow = factor(dow, levels =c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")))
# Plot boxplot manually
Probability <- c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
boxdata <-
  cde4id %>% 
  # filter(id==4669) %>%
  group_by(id, period, dow) %>% 
  summarise(q = quantile(demand, prob = Probability, type=8)) %>%
  mutate(name = rep(paste0("Q", 1:length(Probability)), n_distinct(id, period))) %>% 
  pivot_wider(names_from = name, values_from = q)
boxdata
boxcols <- gghdr:::darken_fill(col = rep("blue", length(Probability)), prob = Probability)
p_box <- boxdata %>% 
  ggplot() +
  geom_rect(aes(ymin=Q1, ymax=Q7, xmin=period-0.2, xmax=period+0.2), color = boxcols[7]) +
  geom_rect(aes(ymin=Q2, ymax=Q6, xmin=period-0.2, xmax=period+0.2), color = boxcols[6]) +
  geom_rect(aes(ymin=Q3, ymax=Q5, xmin=period-0.2, xmax=period+0.2), color = boxcols[5]) +
  geom_rect(aes(ymin=Q4, ymax=Q4, xmin=period-0.5, xmax=period+0.5), color = "black", alpha = 1) +
  # geom_line(aes(x=period, y=Q2), alpha = 0.5) +
  scale_x_continuous(breaks = c(0, 12*(1:4))) + 
  labs(x = "Time of week", y = "Demand (kWh)", color = "Probability") +
  facet_grid(id~dow, scales = "fixed") + 
  scale_color_manual(name = "Probability",
                     values = c("99.0%"=boxcols[7], "95.0%"=boxcols[6], "50.0%"=boxcols[5]), 
                     labels = c("99.0%", "95.0%", "50.0%")) 
p_box
ggsave(paste0("paper/figures/Electricity_2d_plot_hdrbox_7dow_4id_1typical_3anomalous.png"), p_box, width = 10, height = 8, dpi = 500)
save(p_box, file = "data/Electricity_2d_quantileplot_4id.rda")

