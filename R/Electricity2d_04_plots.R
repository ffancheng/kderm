# Same as mlann, to be modified
## Most Up-to-date script to use!
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
ids <- spdemand$id
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

r <- 1
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
method <- ml[1]
load( paste0("data/compareden_electricity_2d_N", N, "_", method, "_radius", radius, "_k", k, "_searchtype", searchtype, "_r", format(r, decimal.mark = "_"), ".rda"))
fisomap <- fisomap
fixden_isomap <- fixden_isomap
p_hdr_isomap <- p_hdr_isomap
p_isomap <- p_isomap
load(paste0("data/metric_", method, "_electricity_2d_radius", radius, "_k", k, "_searchtype", searchtype, ".rda"))
metric_isomap <- metric_isomap
# Y_isomap <- metric_isomap$embedding %>% as.data.frame()
# p_hdr_isomap <- hdrscatterplot_new(Y_isomap[,1], Y_isomap[,2], kde.package = "ks", levels = prob, noutliers = noutliers, label = ids)
# p_isomap <- plot_outlier(x = metric_isomap, gridsize = gridsize, prob = prob, riem.scale = riem.scale, f = fisomap, ell.size = 0)

method <- ml[2]
load( paste0("data/compareden_electricity_2d_N", N, "_", method, "_radius", radius, "_k", k, "_searchtype", searchtype, "_r", format(r, decimal.mark = "_"), ".rda"))
flle <- fisomap
fixden_lle <- fixden_isomap
p_hdr_lle <- p_hdr_isomap
p_lle <- p_isomap
load(paste0("data/metric_", method, "_electricity_2d_radius", radius, "_k", k, "_searchtype", searchtype, ".rda"))
metric_lle <- metric_isomap

method <- ml[3]
load( paste0("data/compareden_electricity_2d_N", N, "_", method, "_radius", radius, "_k", k, "_searchtype", searchtype, "_r", format(r, decimal.mark = "_"), ".rda"))
fle <- fisomap
fixden_le <- fixden_isomap
p_hdr_le <- p_hdr_isomap
p_le <- p_isomap
load(paste0("data/metric_", method, "_electricity_2d_radius", radius, "_k", k, "_searchtype", searchtype, ".rda"))
metric_le <- metric_isomap

method <- ml[4]
load( paste0("data/compareden_electricity_2d_N", N, "_", method, "_radius", radius, "_k", k, "_searchtype", searchtype, "_r", format(r, decimal.mark = "_"), ".rda"))
ftsne <- fisomap
fixden_tsne <- fixden_isomap
p_hdr_tsne <- p_hdr_isomap
p_tsne <- p_isomap
load(paste0("data/metric_", method, "_electricity_2d_radius", radius, "_k", k, "_searchtype", searchtype, ".rda"))
metric_tsne <- metric_isomap

method <- ml[5]
load( paste0("data/compareden_electricity_2d_N", N, "_", method, "_radius", radius, "_k", k, "_searchtype", searchtype, "_r", format(r, decimal.mark = "_"), ".rda"))
fumap <- fisomap
fixden_umap <- fixden_isomap
p_hdr_umap <- p_hdr_isomap
p_umap <- p_isomap
load(paste0("data/metric_", method, "_electricity_2d_radius", radius, "_k", k, "_searchtype", searchtype, ".rda"))
metric_umap <- metric_isomap

cat("Load .rda finished at: "); Sys.time()



###-------------------------------------------------
### Plotting
###-------------------------------------------------

# embedding plot
nolabs <- labs(x = "", y = "")
noytext <-  theme(axis.text.y = element_blank())
p <- (
  ((
       (p_isomap$p + labs(x = "", y = "ISOMAP", title = "DC-KDE") + noytext) |
       (p_hdr_isomap$p + labs(x = "", y = "", title = "KDE") + noytext ) 
  ) & scale_y_continuous(n.breaks = 4)) /
    (( 
         (p_lle$p + labs(x = "", y = "LLE") + noytext) |
         (p_hdr_lle$p + nolabs + noytext) 
    ) & scale_y_continuous(breaks = c(-5, 0, 5))) /
    (( 
         (p_le$p + labs(x = "", y = "Laplacian Eigenmaps") + noytext) |
         (p_hdr_le$p + nolabs + noytext) 
    ) & scale_x_continuous(breaks = c(-.025, 0, .025)) & scale_y_continuous(breaks = c(-.05, 0, .05), limits = c(-.05, .05)) ) /
    (( 
         (p_tsne$p + labs(x = "", y = "t-SNE") + noytext) |
         (p_hdr_tsne$p + nolabs + noytext) 
    ) & scale_y_continuous(n.breaks = 4)) / 
    (( 
         (p_umap$p + labs(x = "", y = "UMAP") + noytext) |
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
ggsave(paste0("paper/figures/Electricity_2d_5levels_outliers_comparison_5ml_r", format(r, decimal.mark = "_"), ".png"), p, width = 8, height = 10, dpi = 300)

###-------------------------------------------------
### Highest density region plot
###-------------------------------------------------
p_hdr_isomap
p_isomap

(
  p_hdr_isomap$p | p_isomap$p
) + 
  plot_layout(guides = "collect") &
  guides(color=guide_legend(nrow=1,byrow=TRUE)) &
  theme(legend.position = 'bottom', 
        axis.title.y = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.title = element_text(hjust = 0.5)) 
s


st <- (((ph1 + labs(title = "Exact NN", x="", y="ISOMAP")) | 
          (pha1 + labs(title = "k-d trees", x="", y="")) |
          (phaa1 + labs(title = "Annoy", x="", y="")))  / 
         ((ph2 + labs(x="", y="LLE")) | 
            (pha2 + labs(x="", y="")) |
            (phaa2 + labs(x="", y=""))) / 
         ((ph3 + labs(x="", y="Laplacian Eigenmaps")) | 
            (pha3 + labs(x="", y="")) |
            (phaa3 + labs(x="", y=""))) / 
         ((ph4 + labs(x="", y="Hessian LLE")) | 
            (pha4 + labs(x="", y="")) |
            (phaa4 + labs(x="", y=""))) /
         ((ph5 + labs(x="", y="t-SNE")) | 
            (pha5 + labs(x="", y="")) |
            (phaa5 + labs(x="", y=""))) /
         ((ph6 + labs(x="", y="UMAP")) | 
            (pha6 + labs(x="", y="")) |
            (phaa6 + labs(x="", y="")))
) + 
  plot_layout(guides = "collect") &
  guides(color = guide_legend(nrow = 1,byrow = TRUE)) &
  theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5, face = "bold")) 
st

# ((ph1 / pha1 / phaa1) | (ph2 / pha2 / phaa2) | (ph3 / pha3 / phaa3) | (ph4 / pha4 / phaa4) | (ph5 / pha5 / phaa5) | (ph6 / pha6 / phaa6)) +
#   plot_layout(guides = "collect") &
#   theme(legend.position = 'bottom') &
#   labs(x = "", y = "")

## hdrXX_compare_1id336tow.png
# ggsave(paste0(folder, "hdr10_compare4ml_kdtreeannoy_1id_6methods", filename, ".png"), s, width=10, height=8)
ggsave(paste0(folder, "hdr10_compare4ml_kdtreeannoy_allids_nt", nt, "_", filename, ".png"), s, width=10, height=6, dpi = 500)
ggsave(paste0(folder, "hdr10_compare4ml_kdtreeannoy_allids_nt", nt, "_", filename, "_transpose.png"), st, width=8, height=11, dpi = 500)


( (ph1 + labs(x="ISO1", y="ISO2")) | (ph2 + labs(x="LLE1", y="LLE2")) )+
  plot_layout(guides = "collect") &
  guides(color=guide_legend(nrow=1,byrow=TRUE)) &
  theme(legend.position = 'bottom') -> hdr_isomap_lle
hdr_isomap_lle
ggsave(paste0(folder, "hdr10_compareisomaplle_kdtreeannoy_nt", nt, "_", filename, ".png"), hdr_isomap_lle, width=10, height=6)

ggsave(paste0("~/Desktop/hdr10_isomap_kdtreeannoy_nt", nt, "_", filename, ".png"), width=8, height=6)




# For single household all tows
###-------------------------------------------------
## plot typical and anomalous tow
###-------------------------------------------------
# The most typical ones are those with Region=="1" in yellow.
# The most anomalous ones are those with Region==">99" in black. 
anomaly_nn <- ph1$data %>%  # View()
  rownames_to_column("tow") %>% 
  rename(ISO1 = Y_isomap...1., ISO2 = Y_isomap...2.) %>% 
  mutate(tow = as.numeric(tow), 
         dow = ceiling(tow/48),
         tod = (tow - 48*floor(tow/48))/2,
         method = "nn") %>% 
  filter(Region %in% c(">99", "1"))

anomaly_kdtree <- pha1$data %>%  # View()
  rownames_to_column("tow") %>% 
  rename(ISO1 = Y_annisomap...1., ISO2 = Y_annisomap...2.) %>%
  mutate(tow = as.numeric(tow), 
         dow = ceiling(tow/48),
         tod = (tow - 48*floor(tow/48))/2,
         method = "kdtree") %>% 
  filter(Region %in% c(">99", "1"))

anomaly_annoy <- phaa1$data %>%  # View()
  rownames_to_column("tow") %>% 
  rename(ISO1 = Y_annisomap_annoy...1., ISO2 = Y_annisomap_annoy...2.) %>%
  mutate(tow = as.numeric(tow), 
         dow = ceiling(tow/48),
         tod = (tow - 48*floor(tow/48))/2,
         method = "annoy") %>% 
  filter(Region %in% c(">99", "1"))

bind_rows(anomaly_nn, anomaly_kdtree, anomaly_annoy)
# tow       ISO1        ISO2  Region dow  tod method
# 1   24 -0.1695829 0.2909404      1   1 12.0     nn
# 2   77 -0.1179890 0.3179422      1   2 14.5     nn
# 3   81 -0.1458833 0.2921770      1   2 16.5     nn
# 4  310 -0.1476710 0.3262820      1   7 11.0     nn
# 5   38  0.9525183 0.4212948    >99   1 19.0     nn
# 6  134  1.0237164 0.3636394    >99   3 19.0     nn
# 7  182  0.9442243 0.4780956    >99   4 19.0     nn
# 8  327  1.0776950 0.2679600    >99   7 19.5     nn
# 9   24 -0.1865536 0.2911935      1   1 12.0 kdtree
# 10  77 -0.1344318 0.3245840      1   2 14.5 kdtree
# 11  81 -0.1611130 0.2947367      1   2 16.5 kdtree
# 12 310 -0.1664895 0.3296634      1   7 11.0 kdtree
# 13 134  1.0078732 0.3874026    >99   3 19.0 kdtree
# 14 135  1.1850625 0.1991039    >99   3 19.5 kdtree
# 15 327  1.0640752 0.3013102    >99   7 19.5 kdtree
# 16 328  1.1531247 0.1952406    >99   7 20.0 kdtree
# 17  24 -0.1585714 0.2859572      1   1 12.0  annoy
# 18  30 -0.1534337 0.2466920      1   1 15.0  annoy
# 19 222 -0.2069617 0.2873975      1   5 15.0  annoy
# 20 310 -0.1330758 0.3084639      1   7 11.0  annoy
# 21 134  1.0927996 0.3916109    >99   3 19.0  annoy
# 22 182  0.9507067 0.4667087    >99   4 19.0  annoy
# 23 326  0.9436023 0.4970019    >99   7 19.0  annoy
# 24 327  1.1539380 0.2926905    >99   7 19.5  annoy

# Typical: Wed 19:00, tow 134, 327
# Anomalys: Mon 12:00, tow 24, 310

# Plot the electricity demand distribution for two time periods, index 134 and 24
load("data/DT.rda")
head(DT)
DT_1id <- DT[id == "1003",]
# tows <- c(134, 24, 327, 310) # single household, compare tow
tows <- c(32, 119, 37, 325) # single household, count_ratio
DT2tow <- DT_1id[tow %in% tows,][, -"id"]
DT2tow
ptow <- DT2tow %>% 
  as_tibble() %>% 
  mutate(tow = as.factor(tow),
         typical = ifelse(tow %in% c(32, 119), TRUE, FALSE)) %>% # c(24, 310)
  ggplot(aes(x=demand, group=tow, fill=typical)) + 
  geom_histogram(
    aes(y=..density..), # Histogram with density instead of count on y-axis
    # binwidth=.1,
    bins = 50,
    colour="grey40", alpha=.5, # position="density"
  ) + 
  # geom_density(alpha=.3) +  # Overlay with transparent density plot, removed since fitted using only positive values, not accurate
  scale_fill_manual(values = c("black", "orange")) +
  facet_wrap(~tow, nrow = 2,
             # labeller = as_labeller(c("24"="Monday 12pm", "134"="Wednesday 7pm", "310"="Sunday 11am", "327"="Sunday 11pm"))
             labeller = as_labeller(c("32"="Monday 4 pm", "37"="Monday 6:30 pm", "119"="Wednesday 11:30 am", "325"="Sunday 6:30 pm"))) +
  labs(x = "Electricity demand (kWh)", y = "Density", fill = "Time of week") + 
  theme(legend.position = "None")
ptow
ggsave(paste0(folder, "electricity_compare2tow_1id336tow.png"), ptow, width=10, height=6, dpi = 500)


# # Plot pmf
# load("data/spdemand_3639id_336tow_201length.rda")
# spd2tow <- spdemand[id == "1003" & (tow %in% tows),][, -"id"]
# spd2tow
# pdtow <- spd2tow %>%
#   gather(key = "p", value = "prob", -tow) %>%
#   mutate(p = as.numeric(str_remove(p, "p"))) %>%
#   ggplot(aes(x = p, y = prob, group = tow, col = factor(tow))) +
#   geom_line(aes(linetype = factor(tow))) +
#   # facet_grid(tow ~ .) +
#   labs(x = "", y = "Probability", col = "Time of week", linetype = "Time of week")
# pdtow


# ###-------------------------------------------------
# ### plot typical households (run once)
# ###-------------------------------------------------
# # Plot the electricity demand for two sample IDs, 1003 and 1539
# # load("data/DT.rda")
# # head(DT)
# # DT2id <- DT[id == 1003 | id == 1539,]
# # save(DT2id, file = "data/DT2id.eda")
# load("data/DT2id.rda")
# DT2id
# summary(DT2id)
# p <- DT2id %>%
#   ggplot(aes(x = day, y = demand, group = id)) +
#   geom_line() +
#   facet_grid(id ~ .) +
#   labs(x = "Days", y = "Demand (kWh)")
# ggsave(folder, "smartmeter.png", p, width = 8, height = 6)
# plotly::ggplotly(p, width = 600, height = 400)



# For all households
# id index 451 and 2425 to compare
###-------------------------------------------------
### plot typical households (run once)
###-------------------------------------------------
# Find typical and anomalous households
# The most typical ones are those with Region=="1" in yellow.
# The most anomalous ones are those with Region==">99" in black. 
anomaly_nn <- ph1$data %>%  # View()
  rownames_to_column("index") %>%
  rename(ISO1 = Y_isomap...1., ISO2 = Y_isomap...2.) %>% 
  mutate(
    index = as.numeric(index),
    ID = ids[index],
    method = "nn") %>%
  filter(Region %in% c(">99", "1"))

anomaly_kdtree <- pha1$data %>%  # View()
  rownames_to_column("index") %>%
  rename(ISO1 = Y_annisomap...1., ISO2 = Y_annisomap...2.) %>% 
  mutate(
    index = as.numeric(index),
    ID = ids[index],
    method = "nn") %>%
  filter(Region %in% c(">99", "1"))

anomaly_annoy <- phaa1$data %>%  # View()
  rownames_to_column("index") %>%
  rename(ISO1 = Y_annisomap_annoy...1., ISO2 = Y_annisomap_annoy...2.) %>% 
  mutate(
    index = as.numeric(index),
    ID = ids[index],
    method = "nn") %>%
  filter(Region %in% c(">99", "1"))

bind_rows(anomaly_nn, anomaly_kdtree, anomaly_annoy)
ph1[["plot_env"]][["label"]] %>% head(n=10) # set noutliers = 3639
ph1[["plot_env"]][["label"]] %>% tail(n=10) 
# 2838, 7324, 5136 are anomalies from Isomap, UMAP
# pick 6762 as typical from Isomap
which(ids %in% c(6762, 2838, 7324)) # index of IDs

# DT2id <- DT[id %in% c(6762, 2838, 7324),] # 1321 typical
DT2id <- DT[id %in% c(7049, 5136, 1472),]
save(DT2id, file = "data/electricityplot_DTcompare3id.rda")

# Plot the electricity demand for two sample IDs, index 451 and 2425
library(tidyverse)
# load("data/DT.rda")
# head(DT)
# ids <- DT %>% dtplyr::lazy_dt() %>% pull(id) %>% unique() %>% sort()
# # DT2id <- DT[id %in% c(ids[c(485, 1280)], 1003),] # 4273 typical
# id2 <- ids[c(2057, 1280, 169)]
# DT2id <- DT[id %in% id2,] # 1321 typical
# save(DT2id, file = "data/electricityplot_DTcompare3id.rda")

load("data/electricityplot_DTcompare3id.rda")
DT2id
summary(DT2id)
p <- DT2id %>%
  ggplot(aes(x = day, y = demand, group = id)) +
  geom_line() +
  facet_grid(factor(id, levels = c(1472, 7049, 5136)) ~ .) +
  labs(x = "Days", y = "Demand (kWh)")
p
ggsave(paste0(folder, "electricity_compare3id.png"), p, width = 8, height = 6, dpi = 500)
save(p, file = "data/electricityplot_compare3id.rda")

# Plot the electricity demand distribution for two anomalous households, index 481 and 1280,  and 1 typical household, index 169
# load("data/DT.rda")
# head(DT)
# DT_1id <- DT[id == "1003",]
# tows <- c(134, 24, 327, 310) # single household, compare tow
# DT2tow <- DT_1id[tow %in% tows,][, -"id"]
DT2id
p3ids <- DT2id %>% 
  as_tibble() %>% 
  mutate(id = as.factor(id),
         typical = ifelse(id == unique(DT2id$id)[1], TRUE, FALSE)) %>% 
  ggplot(aes(x=demand, group=id, fill=typical)) + 
  geom_histogram(
    aes(y=..density..), # Histogram with density instead of count on y-axis
    binwidth=.1,
    colour="grey40", alpha=.2, # position="density"
  ) + 
  geom_density(alpha=.3) +  # Overlay with transparent density plot
  scale_fill_manual(values = c("black", "orange")) +
  facet_wrap(~id, ncol = 1, scales = "free_y", strip.position="right") + # ,labeller = as_labeller(c("4669"="ID 4669", "3243"="ID 3243", "1321"="ID 1321"))
  labs(x = "Demand (kWh)", y = "Density") + 
  theme(legend.position = "None")
p3ids

# wrap_plots(p, p3ids)
pcomp <- cowplot::plot_grid(
  p, p3ids,
  align = "h", axis = "tb",
  nrow = 1, rel_widths = c(2, 1)
)
pcomp
ggsave(paste0(folder, "electricity_compare3id_isomap.png"), pcomp, width = 8, height = 6)
# plotly::ggplotly(p, width = 600, height = 400)





## KDE and DC-KDE raw time series plot
load("data/electricityplot_compare3id.rda")
p_dc <- p
load("data/electricityplot_DTcompare3id.rda")
((p_dc + ggtitle("DC-KDE")) +
    (p + ggtitle("KDE"))) + 
  plot_layout(guides = "collect") &
  theme(#legend.position = 'bottom', 
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12)) 
ggsave("paper/figures/dckde_2densities_raw_electricity.png")



###-------------------------------------------------
### plot 3 typical households distributions
###-------------------------------------------------
# # Or: plot the distribution of 336 tows for each household over 1.5 years
# # cde branch of gghdr github package, +geom_hdr_boxplot(aes(x=tow, y=demand))
# # remotes::install_github("ropenscilabs/gghdr@cde")
# remotes::install_github("ffancheng/gghdr") # add nxmargin argument for number of boxplots
# library(ggplot2)
# library(gghdr)
# ggplot(faithful, aes(y = eruptions)) +
#   geom_hdr_boxplot()
# ggplot(faithful, aes(x = waiting, y = eruptions)) +
#   geom_hdr_boxplot(nxmargin=20)
# ggplot(faithful) +
#   geom_point(aes(x = eruptions, y = waiting)) + 
#   geom_hdr_rug(aes(x = eruptions), prob = c(0.99, 0.5), fill = "blue")
# ggplot(data = faithful, aes(x = waiting, y=eruptions)) +
#   geom_point(aes(colour = hdr_bin(x = waiting, y = eruptions))) +
#   scale_colour_viridis_d(direction = -1) 
# 
# # hdrcde examples
# library(hdrcde)
# hdr.boxplot(faithful$eruptions)
# faithful.cde <- cde(faithful$waiting, faithful$eruptions,
#                     x.name="Waiting time", y.name="Duration time")
# plot(faithful.cde)
# plot(faithful.cde,xlab="Waiting time",ylab="Duration time",plot.fn="hdr")
# hdrscatterplot(faithful$waiting, faithful$eruptions)
# 
# 
# # hdrcde for 3 ids
# library(tidyverse)
# library(hdrcde)
# load("data/electricityplot_DTcompare3id.rda")
# cde2id <-
#   DT2id %>%
#     as_tibble() %>%
#     mutate(tow = as.integer(tow))
# 
# 
# # id.cde0 <- cde(cde2id$tow, cde2id$demand, x.margin = 1:336,              x.name = "Time of week period", y.name = "Demand")
# # Or equivalent nxmargin=336
# id.cde <- cde(cde2id$tow, cde2id$demand, nxmargin = 336,
#               x.name = "Time of week period", y.name = "Demand")
# plot(id.cde)
# plot(id.cde, plot.fn="hdr")
# 
# par(mfrow=c(3,1))
# for (i in 1:3) {
#   cdeplot <- filter(as_tibble(cde2id), id == unique(cde2id$id)[i])
#   cde(cdeplot$tow, cdeplot$demand, x.margin = 1:336) %>%
#     plot(xlab="Time of week period", ylab="Demand", plot.fn="hdr")
# }


# modified gghdr for 3 ids
# remotes::install_github("ffancheng/gghdr") # add nxmargin argument for number of boxplots
library(ggplot2)
library(gghdr)
load("data/electricityplot_DTcompare3id.rda")
cde2id <- 
  DT2id %>% 
  as_tibble() %>% 
  mutate(tow = as.integer(tow),
         dow = str_replace_all(dow, c("1"="Monday", "2"="Tuesday", "3" = "Wednesday", "4" = "Thursday", "5" = "Friday", "6" = "Saturday", "7" = "Sunday")), 
         dow = factor(dow, levels =c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")))
# # cde2id <- cde2id[c(1:100, 25726+(1:100)),]
# # # plot against 336 time of the week
# # p_idtow <- cde2id %>%
# #   ggplot(aes(x=tow, y=demand)) +
# #   geom_hdr_boxplot(nxmargin=336, fill = "blue") +
# #   facet_wrap(~id, ncol = 1) +
# #   theme(legend.position = 'bottom') +
# #   labs(x = "Time of week", y = "Demand (kWh)", fill = "Probability", colour = "Probability")
# # ggsave(folder, "electricity_gghdr_3id_336tow.png", width = 8, height = 6)
# # plot against 48 period of day and facet by 7 day of the week
# p <- cde2id %>%
#   ggplot(aes(x=period, y=demand)) +
#   geom_hdr_boxplot(nxmargin=48, fill = "blue")
# p_iddow <- p +
#   scale_x_continuous(breaks = c(0, 12*(1:4))) +
#   facet_grid(id~dow, 
#              # labeller = as_labeller(c("1"="Monday", "2"="Tuesday", "3" = "Wednesday", "4" = "Thursday", "5" = "Friday", "6" = "Saturday", "7" = "Sunday"))
#              ) +
#   # theme(strip.text.x = element_text(size = 5)) + 
#   labs(x = "Time of week", y = "Demand (kWh)", fill = "Probability", colour = "Probability") + # Legend title not changing?
#   # scale_color_grey(name = "Prob") +
#   # guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
#   theme(legend.position = 'bottom')
# p_iddow
# ggsave(folder, "electricity_gghdr_3id_7dow.png", p_iddow, width = 12, height = 8)
# 



# Plot boxplot manually
Probability <- c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
boxdata <-
  cde2id %>% 
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
  facet_grid(factor(id, levels = c(1472, 7049, 5136))~dow, scales = "fixed") + 
  scale_color_manual(name = "Probability",
                     values = c("99.0%"=boxcols[7], "95.0%"=boxcols[6], "50.0%"=boxcols[5]), 
                     labels = c("99.0%", "95.0%", "50.0%")) 
p_box
ggsave(paste0(folder, "electricity_hdrbox_3id_7dow.png"), p_box, width = 12, height = 8, dpi = 500)
save(p_box, file = "data/electricity_quantileplot_3ids.rda")





###-------------------------------------------------
## specificity and sensitivity
###-------------------------------------------------
# xtab_set <- function(A, B){
#   both <- union(A,B)
#   KNN <- both %in% A %>% factor(levels=c(TRUE, FALSE))
#   ANN <- both %in% B %>% factor(levels=c(TRUE, FALSE))
#   return(table(KNN, ANN))
# }
# # set.seed(1)
# # A <- sample(letters[1:20],10,replace=TRUE)
# # B <- sample(letters[1:20],10,replace=TRUE)
# # xtab_set(A,B)
# # #        inB
# # # inA     FALSE TRUE
# # #   FALSE     0    5
# # #   TRUE      6    3
# 
# ht1 <- xtab_set(A=ph1$plot_env$outliers, B=pha1$plot_env$outliers)
# ht2 <- xtab_set(A=ph2$plot_env$outliers, B=pha2$plot_env$outliers)
# ht3 <- xtab_set(A=ph3$plot_env$outliers, B=pha3$plot_env$outliers)
# ht4 <- xtab_set(A=ph4$plot_env$outliers, B=pha4$plot_env$outliers)
# 
# library(caret)
# specificity(ht1) # Specificity is always 0, because the outlier union appears at least in one set
# sensitivity(ht1)
# sensitivity(ht2)
# sensitivity(ht3)
# sensitivity(ht4)
# 
# # For ID 1003
# # > sensitivity(ht1)
# # [1] 0.45
# # > sensitivity(ht2)
# # [1] 0.4
# # > sensitivity(ht3)
# # [1] 0.6
# # > sensitivity(ht4)
# # [1] 0.55


