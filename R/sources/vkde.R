# input x is a multivariate data/matrix
# takes care of the computation when bandwidth is full matrix instead of fixed vector or diagonal matrix
# lims is a 2-column matrix of ranges(xmin, xmax) for each column in x
vkde <- function(x, h, gridsize = 25, xmin = apply(x, 2, min), xmax = apply(x, 2, max), eval.points){

  n <- nrow(x)
  d <- ncol(x)
  if(any(!is.finite(x)))
    stop("missing or infinite values in the data are not allowed")
  if(any(!is.finite(xmin) | !is.finite(xmax)))
    stop("only finite values are allowed in the minimum/maximum values for grids")
  if(is_scalar_atomic(gridsize)) gridsize <- rep(gridsize, d) # vector of number of grid points

  if(is.vector(h)) return(ks::kde(x, h = h, gridsize = gridsize, xmin = xmin, xmax = xmax, eval.points = eval.points)) # fixed diagonal bandwidth # return(MASS::kde2d(x, y, h, n, lims)) # only works for 2d

  if(missing(eval.points)) {
    # grid points for multidimension, linear method
    gx <- numeric(0)
    for (i in 1:d) {
      gx <- c(gx, list(seq(xmin[i], xmax[i], length.out = gridsize[i])))
    }
    g <- expand.grid(gx)
    
  } else {
    g <- eval.points
  }
  
  p <- prod(gridsize)
  # bandwidth is a d*d*nx array
  # z <- array(NA, gridsize)
  z <- NULL
  
  for (k in 1:n) {
    hk <- h[,,k]
    z <-  abind::abind(z, array(mvtnorm::dmvnorm(x = g, mean = x[k,], sigma = hk), dim = gridsize), along = d + 1)  # stack array of dimension (gridsize*gridsize) with abind
  }
  z <- rowMeans(z, dims = 2, na.rm = TRUE)
  
  ## TODO: optimize hk as bandwidth pointwise
  # AMISE
  
  f <- c(x = gx, list(z = z))
  if(d == 2) names(f)[1:2] <- c("x", "y")
  return(f)
  
}
