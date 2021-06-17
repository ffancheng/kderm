# Variable kernel density estimate
# \hat{f}(x) = n^{-1} \sum_i K_{H_i} (x - X_i)

# calculate the density estimate for each point in emb_isomap with each H matrix, H_isomap[i,,]

ker.density <- function(x, h) {
  x <- sort(x)
  n <- length(x)
  s <- 0
  t <- 0
  y <- 0
  for (i in 2:n)
    s[i] <- 0
  for (i in 1:n) {
    for (j in 1:n)
      s[i] <- s[i] + exp(-((x[i] - x[j]) ^ 2) / (2 * h * h)) # h^2?    
      t[i] <- s[i]
  }
  for (i in 1:n)
    y[i] <- t[i] / (n * h * sqrt(2 * pi))
  z <- complex(re = x, im = y)
  hist(x, freq = FALSE)
  lines(z)
}

# # univariate
# dfn1 <- function(x) {
#   0.5 * dnorm(x, 3, 1) + 0.5 * dnorm(x, -3, 1)
# }
# par(mfrow = c(2, 2))
# curve(dfn1(x), from = -6, to = 6)
# data <- c(rnorm(200, 3, 1), rnorm(200, -3, 1))
# plot(density(data,bw=8))
# plot(density(data,bw=0.8))
# plot(density(data,bw=0.08))
# 
# par(mfrow = c(1,1))
# ker.density(data, 0.8)


# 2-d kde
mkde <- function (x, h) {
  # Data is a matrix of n*d
  # H is an array of dimension (d,d,n)
  # start <- Sys.time()

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
      s[i] <- s[i] + det(h[,,i])^(-1/2) * exp(- 1/2 * t(x[i,] - x[j,]) %*% solve(h[,,i]) %*% (x[i,] - x[j,]))
    }
    y[i] <- s[i] / (n * 2 * pi)
  }
  # print(Sys.time() - start)

  return(y)
}
# x = metric_isomap$embedding
# f <- mkde(x = metric_isomap$embedding, h = riem_isomap)
# # is.na(f) <- sapply(f, is.infinite)
# # matrixcalc::is.positive.definite(riem_isomap[,,1])


# multivariate variable density estimate
# https://www.mathworks.com/matlabcentral/fileexchange/68762-2d-kernal-density-estimation
# https://github.com/cran/MASS/blob/master/man/kde2d.Rd

# Refer to both the R and MATLAB code about multivariate kernel density estimate
# x_i is each original data point
# x is the constructed grid points to calculate the density estimate
# 
vkde2d <- function(x, y, h, n = 25, lims = c(range(x), range(y)) ){
  
  nx <- length(x)
  if(length(y) != nx)
    stop("data vectors must be the same length")
  if(any(!is.finite(x)) || any(!is.finite(y)))
    stop("missing or infinite values in the data are not allowed")
  if(any(!is.finite(lims)))
    stop("only finite values are allowed in 'lims'")

  # grid points
  n <- rep(n, length.out = 2L)  # number of grids for all dimensions
  gx <- seq.int(lims[1L], lims[2L], length.out = n[1L])
  gy <- seq.int(lims[3L], lims[4L], length.out = n[2L])
  
  # bandwidth matrix is a n*2 matrix 
  h <- if (missing(h)) t(replicate(nx, c(bandwidth.nrd(x), bandwidth.nrd(y)))) else h
        # else rep(h, length.out = 2L) # for fixed h input in kde2d
          
  if (any(h <= 0))
    stop("bandwidths must be strictly positive")
  
  # z <- matrix(0, n[1L], n[2L])
  # 
  # for(i in 1:n[1L]) {
  #   for (j in 1:n[2L]) {
  #     
  #     s = 0
  #     for (k in 1:nx) {
  #       s <- s + 1 / (2 * pi) * (-1/(2*h[k,1]*h[k,2])) * exp(- 1/2 * t(gx[i] - x[k]) * 1/(h[k, 1]*h[k, 2]) * (gy[j] - y[k]))
  #     }
  #     z[i, j] <- s / nx
  #     
  #   }
  # }
  
  h <- h/4                            # for S's bandwidth scale
  ax <- outer(gx, x, "-" )/h[1L]
  ay <- outer(gy, y, "-" )/h[2L]
  z <- tcrossprod( matrix(dnorm(ax), , nx)/h[1L], matrix(dnorm(ay), , nx)/h[2L] ) /
    nx
  # z <- tcrossprod(matrix(dnorm(ax), , nx), matrix(dnorm(ay), , nx)) /
  #   (nx * h[1L] * h[2L])
  
  return(list(x = gx, y = gy, z = z))
}

# attach(geyser)
# plot(duration, waiting, xlim = c(0.5,6), ylim = c(40,100))
# f1 <- kde2d(x = duration, y = waiting, n = 50, lims = c(0.5, 6, 40, 100))
# image(f1, zlim = c(0, 0.05))
# f2 <- vkde2d(x = duration, y = waiting, n = 50, lims = c(0.5, 6, 40, 100))
# image(f2, zlim = c(0, 0.05))
