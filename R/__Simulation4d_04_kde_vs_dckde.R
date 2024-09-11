## -----------------------------------------------------------------------------
## This script is for manifold learning and kernel density estimation (KDE and DC-KDE)
## -----------------------------------------------------------------------------
# Modify `method` as different manifold learning methods in embed()
# `radius=10`
# Tune `r`, the bandwidth parameter, for DC-KDE
# One ggplot for rank plots, and two .rda files for ML and density estimates
# This script is modified for HPC jobs for four ML methods
rm(list = ls())
library(tidyverse)
library(dimRed)
library(hdrcde)
library(ggforce)
library(ks)
library(patchwork)
Jmisc::sourceAll(here::here("R/sources"))
set.seed(1234)
# r <- 1
# k <- 500
load(paste0("data/simdata_100d_4dmanifold_N10000_trueden_k", k, ".rda"))


# ----KDE with ks::kde on high-dimensional data instead of embedding------------
library(weird)
library(tictoc)

# From Rob
n <- 100
d <- 3
y <- matrix(rnorm(n*d), ncol = d)
H <- kde_bandwidth(y)
fhat <- kde(y, H)


# set H=diag(d)
# X: true manifold, train: 100d data for dimRed, N, p, den: density of 5d semishepre, trueden: true density of 4d manifold

# KDE on true manifold X
dim(X) # 10000*4
dim(train) # 10000*100
tic()
H <- kde_bandwidth(X)
toc() # 0.141s
fhat <- kde(X, H)
summary(fhat$estimate)
summary(trueden)

# set H as Identity matrix for 100d data
H1 <- diag(ncol(train)) #+ diag(1e-8, ncol(train))
tic()
fhat1 <- kde(x = train, H = H1, binned = FALSE, eval.points = train)
toc() # 20min
summary(fhat1$estimate)
summary(trueden)



# ----Compare with true density and kde, dckde----------------------------------

