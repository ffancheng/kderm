library(tidyverse)
library(dimRed)
library(dtplyr)
library(reticulate)
library(viridis)
library(ks)

load("data/spdemand_3639id336tow.rda")
nid <- 1
ntow <- 336

train <- spdemand %>%
  lazy_dt() %>%
  filter(tow <= ntow,
         # id <= sort(unique(spdemand[,id]))[nid],
         id == 1003) %>%
  dplyr::select(-id, -tow) %>%
  as.data.table()
N <- nrow(train)

# Parameters fixed
D <- 2
K <- 20

# RColorBrewer::display.brewer.all()
# cols <- RColorBrewer::brewer.pal(6, "Blues")


# Check the embedding plot using R and Python
# Isomap looks similar

# dimRed embedding plot
png("figures/dimred_embedding.png",
    width = 800,
    height = 600)
par(mfrow = c(2, 2))
for (i in 1:3) {
  method <- c("Isomap", "LLE", "LaplacianEigenmaps")[i]
  e <-
    dimRed::embed(
      train,
      .method = method,
      knn = K,
      ndim = D,
      .mute = c("output")
    )
  # plot(e, type = "2vars")
  plot(e@data@data, col = viridis::viridis(24))
}
dev.off()

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
# hat(f)(x) = n^(-1) sum_i K_H (x - X_i).

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
      s[i] = s[i] + exp(-((x[i] - x[j]) ^ 2) / (2 * h * h))
    t[i] = s[i]
  }
  for (i in 1:n)
    y[i] = t[i] / (n * h * sqrt(2 * pi))
  z = complex(re = x, im = y)
  hist(x, freq = FALSE)
  lines(z)
}
ker.density(data, 0.8)


# 2-d kde
kde2d_mod <- function (data, h) {
  # Data is a matrix
  print(Sys.time()) #for timing
  
  nx <- dim(data)[1]
  if (dim(data)[2] != 2) 
    stop("data vectors have only one variable")
  if (any(!is.finite(data))) 
    stop("missing or infinite values in the data are not allowed")
  if (any(!is.finite(lims))) 
    stop("only finite values are allowed in 'lims'")
  
  for (i in 2:n)
    s[i] = 0
  for (i in 1:n) {
    for (j in 1:n)
      s[i] = s[i] + exp(-((x[i] - x[j]) ^ 2) %*% solve(2 * h %*% h))
    t[i] = s[i]
  }
  for (i in 1:n)
    y[i] = t[i] / (n * h * sqrt(2 * pi))
  z = complex(re = x, im = y)
  hist(x, freq = FALSE)
  lines(z)
  }
  
  print(Sys.time())
  #Un-transpose the final data.
  z<-t(matrix(z,n,n))
  dim(z)<-c(n^2,1)
  z<-as.vector(z)
  return(z)
}

kde2d_mod(data = as.matrix(emb_isomap[1, ]), h = H_isomap[1,,])
          

# kernel estimated in grid
grid<- function(n,lims) {
  num <- rep(n, length.out = 2L)
  gx <- seq.int(lims[1L], lims[2L], length.out = num[1L])
  gy <- seq.int(lims[3L], lims[4L], length.out = num[2L])
  
  v1=rep(gy,length(gx))
  v2=rep(gx,length(gy))
  v1<-matrix(v1, nrow=length(gy), ncol=length(gx))
  v2<-t(matrix(v2, nrow=length(gx), ncol=length(gy)))
  grid_out<-c(unlist(v1),unlist(v2))
  
  grid_out<-aperm(array(grid_out,dim=c(n,n,2)),c(3,2,1) ) #reshape
  grid_out<-unlist(as.list(grid_out))
  dim(grid_out)<-c(2,n^2)
  grid_out<-t(grid_out)
  return(grid_out)
}

# plotting
kde2d_mod_plot<-function(kde2d_mod_output,n,lims){
  num <- rep(n, length.out = 2L)
  gx <- seq.int(lims[1L], lims[2L], length.out = num[1L])
  gy <- seq.int(lims[3L], lims[4L], length.out = num[2L])
  
  v1=rep(gy,length(gx))
  v2=rep(gx,length(gy))
  v1<-matrix(v1, nrow=length(gy), ncol=length(gx))
  v2<-t(matrix(v2, nrow=length(gx), ncol=length(gy)))
  
  image.plot(v1,v2,matrix(kde2d_mod_output,n,n))
  map('world', fill = FALSE,add=TRUE)
}


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
?ks::kde
?ks::vkde

i = 1
for (i in 1:N) {
  f_i <- ks:::kde.balloon.2d(x = as.matrix(emb_isomap[i, ]),
                             H = H_isomap[i, , ])
  plot(f_i, display = "persp")
  
}


data(worldbank)
wb <- as.matrix(na.omit(worldbank[, 4:5]))
xmin <- c(-70, -25)
xmax <- c(25, 70)
fhat <- kde(x = wb, xmin = xmin, xmax = xmax)
fhat.sp <- kde.sp(x = wb, xmin = xmin, xmax = xmax)
plot(
  fhat,
  display = "persp",
  box = TRUE,
  phi = 20,
  zlim = c(0, max(fhat.sp$estimate))
)
plot(
  fhat.sp,
  display = "persp",
  box = TRUE,
  phi = 20,
  zlim = c(0, max(fhat.sp$estimate))
)
if (interactive()) {
  fhat.ball <- kde.balloon(x = wb, xmin = xmin, xmax = xmax)
  plot(
    fhat.ball,
    display = "persp",
    box = TRUE,
    phi = 20,
    zlim = c(0, max(fhat.sp$estimate))
  )
}
