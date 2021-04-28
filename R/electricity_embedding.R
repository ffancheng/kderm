library(tidyverse)
library(dimRed)
library(dtplyr)
library(reticulate)

load("data/spdemand_3639id336tow.rda")
nid <- 1
ntow <- 336

train <- spdemand %>% 
  lazy_dt() %>% 
  filter(tow <= ntow,
    # id <= sort(unique(spdemand[,id]))[nid],
    id == 1003
         ) %>% 
  dplyr::select(-id, -tow) %>% 
  as.data.table() 
N <- nrow(train)

# Parameters fixed
D <- 2
K <- 20

# RColorBrewer::display.brewer.all()
cols <- RColorBrewer::brewer.pal(6, "Blues")
# Check the embedding plot using R and Python
# Isomap looks similar

# dimRed embedding plot
par(mfrow=c(2,2))
for(i in 1:3){
  method <- c("Isomap", "LLE", "LaplacianEigenmaps")[i]
  e <- dimRed::embed(train, .method = method, knn = K, ndim = D, .mute = c("output"))
  # plot(e, type = "2vars")
  plot(e@data@data, col = viridis::viridis(24))
}

# megaman embedding plot
emb_isomap <- feather::read_feather("data/embedding_isomap_1id336tow.feather")
emb_lle <- feather::read_feather("data/embedding_lle_1id336tow.feather")
emb_ltsa <- feather::read_feather("data/embedding_ltsa_1id336tow.feather")
emb_spectral <- feather::read_feather("data/embedding_spectral_1id336tow.feather")
par(mfrow=c(2,2))
plot(emb_isomap$`0`, emb_isomap$`1`, col = viridis::viridis(24))
plot(emb_lle$`0`, emb_lle$`1`, col = viridis::viridis(24))
plot(emb_spectral$`0`, emb_spectral$`1`, col = viridis::viridis(24))
plot(emb_ltsa$`0`, emb_ltsa$`1`, col = viridis::viridis(24))


# H matrix: inverse of Riemmanian matrix, N*d*d array
reticulate::py_config() # manifold_env
np <- import("numpy", convert=TRUE)
H_isomap <- np$load("data/hmatrix_isomap_1id336tow.npy")
H_lle <- np$load("data/hmatrix_lle_1id336tow.npy")
H_ltsa <- np$load("data/hmatrix_ltsa_1id336tow.npy")
H_spectral <- np$load("data/hmatrix_spectral_1id336tow.npy")

H_isomap[1,,] # the matrix for the first data point


# Kernel density estimate



