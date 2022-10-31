rm(list=ls())
library(data.table)
library(dtplyr)
library(tidyverse)
library(ggplot2)
library(dimRed)
library(pryr)
# Jmisc::sourceAll("~/git/kderm/R/sources")
set.seed(1)
nid <- 3639
ntow <- 336
len <- 100
filename <- paste0(nid, "id_", ntow, "tow_", len, "length")

# scen <- 1 # comment on cluster
# scen <- as.numeric(commandArgs()[[6]])

load(paste0('data/half_count_ratio_3639id336tow.rda')) # 3639*(336*100+1)
rownames(spdemand) <- paste0(spdemand$id)
train <- spdemand[, !c("id")]
N <- nrow(train)

distance <- c("euclidean", "manhattan")[2]
print(paste("Distance measure:", distance))
## Parameters
# D <- 7
K <- 100
radius <- as.numeric(commandArgs()[[6]])# 20
treetype <- "kd"
searchtype <- "priority"

## True NN
switch(distance,
       "euclidean" = {truenn <- RANN::nn2(train, query = train, k = K, treetype = treetype, searchtype = searchtype, eps = 0)$nn.idx[, -1]},
       "manhattan" = {truenn <- RANN.L1::nn2(train, query = train, k = K, treetype = treetype, searchtype = searchtype, eps = 0)$nn.idx[, -1]}
)
save(truenn, file = paste0("data/truenn_tvd_", filename, "_radius", radius, "_k", K, "_searchtype", searchtype, ".rda"))
cat("Finished at: "); Sys.time()
