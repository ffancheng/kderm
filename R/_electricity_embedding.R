rm(list = ls())
library(tidyverse)
library(dimRed)
library(dtplyr)
library(data.table)
library(reticulate)
library(viridis)
# library(ks)
library(hdrcde)
library(igraph)
# library(seurat)
library(matrixcalc)
library(intRinsic)
Jmisc::sourceAll("R/sources")
scen <- as.numeric(commandArgs()[[6]])

# load("data/spdemand_3639id336tow.rda")
# nid <- 1
nid <- 3639
ntow <- 336
train <- readRDS(paste0("data/spdemand_", nid, "id", ntow, "tow_train.rds"))

# if(nid == 1) {
#   train <- spdemand %>%
#     lazy_dt() %>%
#     filter(tow <= ntow,
#            # id <= sort(unique(spdemand[,id]))[nid]
#            id == 1003
#     ) %>%
#     dplyr::select(-id, -tow) %>%
#     data.table::as.data.table()
# } else if(nid == 3639) {
#   train <- spdemand %>% 
#     lazy_dt() %>%
#     mutate(id_tow = paste0(id, "_", tow)) %>% 
#     as.data.table() %>%
#     column_to_rownames("id_tow") %>% 
#     select(-c(id, tow))
# }
# # saveRDS(train, file = paste0("data/spdemand_", nid, "id", ntow, "tow_train.rds"))

N <- nrow(train)
# intRinsic::twonn(train, method = c("mle", "linfit", "bayes")[scen]) # 6 for 1ID

# Parameters fixed
x <- train
s <- 2
k <- 100
method <- "annIsomap"
annmethod <- "kdtree"
distance <- c("euclidean", "manhattan")[2]
treetype = "kd"
searchtype = "radius" # change searchtype for radius search based on `radius`, or KNN search based on `k`
radius <- 1
metric_isomap <- metricML(x, s = s, k = k, radius = radius, method = method, annmethod = annmethod, eps = 0, distance = distance, treetype = treetype, searchtype = searchtype) 

# Comparison of Isomap embedding plot
par(mfrow = c(1, 2))
plot(metric_isomap$embedding, col = viridis::viridis(24)) # metricML
# plot(e@data@data, col = viridis::viridis(24)) # embedding from dimRed, same as metricML()
# emb_isomap <-
#   feather::read_feather("data/embedding_isomap_1id336tow.feather")
plot(-emb_isomap$`0`, emb_isomap$`1`, col = viridis::viridis(24)) # embedding from megaman, radius = 0.4

# Use adjacency matrix as input for metricML and dimRed, then the embedding should be close, as well as the Riemannian metric
metric_isomap$adj_matrix # renormalized adjacency matrix
metric_isomap$weighted_graph %>% is.connected()
riem_isomap <- metric_isomap$rmetric
riem_isomap[,,1] %>% 
  # isSymmetric()
  matrixcalc::is.positive.definite() # FALSE

np = import("numpy")
H_isomap <- np$load("data/hmatrix_isomap_1id336tow.npy")
H_isomap[1,,, drop=TRUE] %>% 
  matrixcalc::is.positive.definite() # TRUE

pd <- rep(NA, N)
for (i in 1:N) {
  # pd[i] <- 
  riem_isomap[,,i] %>%  # R
  # H_isomap[i,,, drop=TRUE] %>% # Python
    # isSymmetric()
    # matrixcalc::is.positive.definite() # FALSE
    eigen(only.values = T) %>%
    print()
}
pd %>% sum

determinant(riem_isomap[,,i])
mat <- riem_isomap[,,i]
eigen(mat)

## TODO: plot ellipse to show distortion
emb <- metric_isomap$embedding
i = 1
mat <- riem_isomap[,,i] # not positive definite
center <- c(emb[i,1], emb[i,2])

mat <- matrix(c(2.2, 0.4, 0.4, 2.8), nrow=2) # covariance matrix -> cov(dataMatrix)

library(ellipse)
# plotcorr(mat)
# plotcorr(solve(mat))
plot(ellipse::ellipse(mat, center = center, scale = c(1,1)), type = "l")

car::ellipse(center, shape=mat, radius=1, col="red", lty=2)
# library(mixtools)
mixtools::ellipse(mu = center, sigma = mat, type = "l")

# Ellipse plot centered at each point with the covariance matrix as the Riem matrix
for(i in 1:N){
  mat <- riem_isomap[,,i]
  center <- c(x[i,1], x[i,2])
  add <- ifelse(i == 1, F, T)
  car::ellipse(center, mat, radius = .002, 
               col = blues9[5], asp = 1, lwd = 1, center.pch = 19, center.cex = 0,
               fill = T, fill.alpha = 0.2, add = add, grid = T,
               xlim = c(-.2, .25), ylim = c(-.2, .2))
}


# save.par <- par(ask = interactive())
# Plot the correlation matrix for the mtcars data full model fit
data(mtcars)
fit <- lm(mpg ~ ., mtcars)
plotcorr(summary(fit, correlation = TRUE)$correlation)
# Plot a second figure with numbers in place of the
# ellipses
plotcorr(summary(fit, correlation = TRUE)$correlation, numbers = TRUE)
# Colour the ellipses to emphasize the differences. The color range
# is based on RColorBrewer's Reds and Blues (suggested by Gregor Gorjanc)
corr.mtcars <- cor(mtcars)
ord <- order(corr.mtcars[1,])
xc <- corr.mtcars[ord, ord]
colors <- c("#A50F15","#DE2D26","#FB6A4A","#FCAE91","#FEE5D9","white",
            "#EFF3FF","#BDD7E7","#6BAED6","#3182BD","#08519C")
plotcorr(xc, col=colors[5*xc + 6])
plotcorr(xc, col=colors[5*xc + 6], type = "upper")
plotcorr(xc, col=colors[5*xc + 6], type = "lower", diag = TRUE)
# par(save.par)

# Same as
# e <- dimRed::embed(x,
#                    .method = method,
#                    knn = k,
#                    ndim = s,
#                    annmethod = annmethod,
#                    radius = radius, 
#                    eps = epsilon, 
#                    nt = nt, 
#                    nlinks = nlinks, 
#                    ef.construction = ef.construction, 
#                    distance = distance,
#                    treetype = treetype,
#                    searchtype = searchtype,
#                    .mute = c("output")
#                    
# )
# plot(e@data@data, col = viridis::viridis(24))

# RColorBrewer::display.brewer.all()
# cols <- RColorBrewer::brewer.pal(6, "Blues")

# Check the embedding plot using R and Python
# Isomap looks similar

# # dimRed embedding plot
# png("figures/dimred_embedding.png",
#     width = 800,
#     height = 600)
# par(mfrow = c(2, 2))
# for (i in 1:3) {
  
  i=1
  method <- c("Isomap", "LLE", "LaplacianEigenmaps")[i]
  e <-
    dimRed::embed(
      x,
      .method = method,
      knn = k,
      ndim = s,
      .mute = c("output")
    )
  # plot(e, type = "2vars")
  plot(e@data@data, col = viridis::viridis(24))
# }
# dev.off()

# megaman embedding plot
emb_isomap <-
  feather::read_feather("data/embedding_isomap_1id336tow.feather")
emb_lle <-
  feather::read_feather("data/embedding_lle_1id336tow.feather")
emb_ltsa <-
  feather::read_feather("data/embedding_ltsa_1id336tow.feather")
emb_spectral <-
  feather::read_feather("data/embedding_spectral_1id336tow.feather")
png("figures/megaman_embedding.png",
    width = 800,
    height = 600)
par(mfrow = c(2, 2))
plot(emb_isomap$`0`, emb_isomap$`1`, col = viridis::viridis(24))
plot(emb_lle$`0`, emb_lle$`1`, col = viridis::viridis(24))
plot(emb_spectral$`0`, emb_spectral$`1`, col = viridis::viridis(24))
plot(emb_ltsa$`0`, emb_ltsa$`1`, col = viridis::viridis(24))
dev.off()

# H matrix: inverse of Riemmanian matrix, N*d*d array
reticulate::py_config() # manifold_env
np <- import("numpy", convert = TRUE)
H_isomap <- np$load("data/hmatrix_isomap_1id336tow.npy")
H_lle <- np$load("data/hmatrix_lle_1id336tow.npy")
H_ltsa <- np$load("data/hmatrix_ltsa_1id336tow.npy")
H_spectral <- np$load("data/hmatrix_spectral_1id336tow.npy")

H_isomap[1, , ] # the matrix for the first data point




# Variable kernel density estimate
# \hat{f}(x) = n^{-1} \sum_i K_H (x - X_i)

# calculate the density estimate for each point in emb_isomap with each H matrix, H_isomap[i,,]

# univariate
dfn1 <- function(x) {
  0.5 * dnorm(x, 3, 1) + 0.5 * dnorm(x, -3, 1)
}
par(mfrow = c(2, 2))
curve(dfn1(x), from = -6, to = 6)
data <- c(rnorm(200, 3, 1), rnorm(200, -3, 1))
plot(density(data,bw=8))
plot(density(data,bw=0.8))
plot(density(data,bw=0.08))

ker.density <- function(x, h) {
  x = sort(x)
  n = length(x)
  s = 0
  t = 0
  y = 0
  for (i in 2:n)
    s[i] = 0
  for (i in 1:n) {
    for (j in 1:n)
      s[i] = s[i] + exp(-((x[i] - x[j]) ^ 2) / (2 * h * h)) # h^2?                      
    t[i] = s[i]
  }
  for (i in 1:n)
    y[i] = t[i] / (n * h * sqrt(2 * pi))
  z = complex(re = x, im = y)
  hist(x, freq = FALSE)
  lines(z)
}
par(mfrow = c(1,1))
ker.density(data, 0.8)


# 2-d kde
mkde <- function (x, h) {
  # Data is a matrix of n*d
  # H is an array of dimension (d,d,n)
  start <- Sys.time()
  
  n <- nrow(x)
  if (dim(x)[2] < 2)
    stop("Data matrix has only one variable.")
  if (any(!is.finite(x))) 
    stop("Missing or infinite values in the data are not allowed! ")
  if (!all.equal(nrow(h), ncol(h), dim(x)[2]))
    stop("The bandwidth matrix is of wrong dimension. ")
  
  s <- rep(0, n)
  y <- rep(0, n)
  for (i in 1:n) {
    for (j in 1:n){
      s[i] <- s[i] + abs(det(h[,,i]))^(-1/2) * exp(- 1/2 * t(x[i,] - x[j,]) %*% solve(h[,,i]) %*% (x[i,] - x[j,]))
    }
    y[i] <- s[i] / (n * 2 * pi)
  }
  print(Sys.time() - start) 

  return(y)
}
x = metric_isomap$embedding
f <- mkde(x = metric_isomap$embedding, h = riem_isomap)
# is.na(f) <- sapply(f, is.infinite)
# matrixcalc::is.positive.definite(riem_isomap[,,1])


# Plotting
embed_den <- as_tibble(cbind(x = x[,1], y = x[,2], z = f)) %>%
  drop_na() # %>% 
  # arrange(x, y)

embed_den_list <- list(x = embed_den$x, y = embed_den$y, z = embed_den$z)

embed_den
hdrscatterplot(embed_den$x, embed_den$y)
hdrinfo <- hdr.2d(embed_den$x, embed_den$y, prob = c(50, 95, 99)) 
plot.hdr2d(hdrinfo)
# If density values are specified
hdrinfo1 <- hdr.2d(embed_den$x, embed_den$y, prob = c(50, 95, 99), den = embed_den_list)
plot.hdr2d(hdrinfo)
filled.contour(x = embed_den$x, y = embed_den$y, z = embed_den$z)



# Contour plot after interpolation
library(akima)
# ?interp()
akima.spl <- interp(embed_den$x, embed_den$y, embed_den$z, nx=100, ny=100, linear=FALSE)
p.zmin <- min(embed_den$z,na.rm=TRUE)
p.zmax <- max(embed_den$z,na.rm=TRUE)
breaks <- pretty(c(p.zmin,p.zmax),10)
# colors <- heat.colors(length(breaks)-1)
hdrscatterplot(embed_den$x, embed_den$y)
plot.hdr2d(hdrinfo)
contour(akima.spl, main = "smooth interp(*, linear = FALSE)", col = viridis(length(breaks)-1), levels=breaks, add=TRUE)
points(embed_den, pch = 20)

filled.contour(akima.spl, color.palette = viridis,
               plot.axes = { axis(1); axis(2);
                 title("smooth  interp(*, linear = FALSE)");
                 points(embed_den, pch = 3, col= hcl(c=20, l = 10))})








# library(gghdr)
# ggplot(data = embed_den, aes(x, y)) +
#   geom_point(aes(colour = hdr_bin(x, y, prob = c(50, 95, 99)/100))) +
#   scale_colour_viridis_d(direction = -1) 
# 
#        
# as_tibble(cbind(x, den = f)) %>%
#   ggplot(aes(E1, E2)) + 
#   geom_point() + 
#   # geom_density2d_filled()
#   # stat_density_2d(aes(fill = ..density..), geom = "raster", contour = "white")
#   stat_density_2d(aes(fill = den), geom = "polygon", colour=FALSE) +
#   scale_fill_distiller(palette=4, direction=-1)
  










###############################################################################
## Multivariate kernel density estimate using normal kernels
##
## Parameters
## x - points
## H - bandwidth matrix
## gridsize - number of interval points in grid
## supp - effective support of kernel
## eval.points - compute density estimate at these points (if missing
##            and dim(x) = 2, 3 compute density estimate over grid)
## eval.levels - compute 3-D in 2-D slices taken at these level curves
##
## Returns
## list with first d components with the points that the density
## estimate is evaluated at, and values of the density estimate
##############################################################################
# ?ks::kde
# ?ks::vkde
# 
# i = 1
# for (i in 1:N) {
#   f_i <- ks:::kde.balloon.2d(x = as.matrix(emb_isomap[i, ]),
#                              H = H_isomap[i, , ])
#   plot(f_i, display = "persp")
#   
# }
# 
# 
# data(worldbank)
# wb <- as.matrix(na.omit(worldbank[, 4:5]))
# xmin <- c(-70, -25)
# xmax <- c(25, 70)
# fhat <- kde(x = wb, xmin = xmin, xmax = xmax)
# fhat.sp <- kde.sp(x = wb, xmin = xmin, xmax = xmax)
# plot(
#   fhat,
#   display = "persp",
#   box = TRUE,
#   phi = 20,
#   zlim = c(0, max(fhat.sp$estimate))
# )
# plot(
#   fhat.sp,
#   display = "persp",
#   box = TRUE,
#   phi = 20,
#   zlim = c(0, max(fhat.sp$estimate))
# )
# if (interactive()) {
#   fhat.ball <- kde.balloon(x = wb, xmin = xmin, xmax = xmax)
#   plot(
#     fhat.ball,
#     display = "persp",
#     box = TRUE,
#     phi = 20,
#     zlim = c(0, max(fhat.sp$estimate))
#   )
# }
