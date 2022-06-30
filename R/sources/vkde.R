# input x is a multivariate data/matrix
# r is the global scalar bandwidth, not the plug-in bandwidth MATRIX from kde
# H is the adaptive bandwidth matrix as a N*d*d array
# h is the bandwidth matrix for kde()
# takes care of the computation when bandwidth is full matrix instead of fixed vector or diagonal matrix
# xmin and xmax are long vectors for all column in x
# eval.points by default make grids and generate estimates of grid points; if provided, eg. x, can evaluate density at certain points
# output long list of x, eval.points (grid points or data points), estimates for density, bandwidth matrix, whether estimates are for grids

vkde <- function(x, h = NULL, vh = NULL, r = NULL, gridsize = 20, xmin = apply(x, 2, min), xmax = apply(x, 2, max), eval.points, xextend = 0.15, kde.package = "ks", positive = FALSE, opt.method = c("AMISE", "MEAN", "SCALED"), riem.scale = 1, adj_matrix = NULL, ...){

  n <- nrow(x)
  d <- ncol(x)
  
  if(any(!is.finite(x)))
    stop("Missing or infinite values in the data are not allowed!")
  if(any(!is.finite(xmin) | !is.finite(xmax)))
    stop("Only finite values are allowed in the minimum/maximum values for grids.")
  if(gridsize <= 0) stop("Gridsize must be positive.")
  if(is_scalar_atomic(gridsize)) gridsize <- rep(gridsize[1], d) # vector of number of grid points

  if(is.null(vh) | !is.array(vh)) return(ks::kde(x, h = h, gridsize = gridsize, xmin = xmin, xmax = xmax, eval.points = eval.points, positive = positive, ...)) # fixed diagonal bandwidth # return(MASS::kde2d(x, y, h, n, lims)) # but only works for 2d
  if(is_scalar_atomic(r)) warning("Input r is not a scalar! Default value r=1 is selected."); r <- 1
  
  # Use optimized bandwidth from minimizing AMISE as r in Pelletier's estimator
  if (d == 1 & !positive) 
    h <- ks::hpi(x = x, nstage = 2, binned = ks:::default.bflag(d = d, n = n), deriv.order = 0)
  if (d > 1 & !positive) 
    h <- ks::Hpi(x = x, nstage = 2, binned = ks:::default.bflag(d = d, n = n), deriv.order = 0)
  
  opt.method <- match.arg(opt.method, c("AMISE", "MEAN", "SCALED"), several.ok = FALSE)
  xr <- apply(x, 2, function(x) diff(range(x, na.rm = TRUE)))
  xmin <- xmin - xr * xextend
  xmax <- xmax + xr * xextend

  # bandwidth is a d*d*nx array
  ## TODO: optimize hi as bandwidth pointwise
  if(opt.method == "MEAN"){
    # Option 1: scale hi with mean_i(|hi|) / det(hi)
    vhidet <- apply(vh, 3, det)
    vh <- sweep(vh, 3, mean(hidet) / hidet, "*")
  } else
  
  if(opt.method == "AMISE"){
  # Option 2: scale hi as mean(det(hi)) * const = det(hopt) 
  # # Use optimized bandwidth from minimizing AMISE for scaling
  # if (d == 1 & !positive) 
  #   H <- ks::hpi(x = x, nstage = 2, binned = ks:::default.bflag(d = d, n = n), deriv.order = 0)
  # if (d > 1 & !positive) 
  #   H <- ks::Hpi(x = x, nstage = 2, binned = ks:::default.bflag(d = d, n = n), deriv.order = 0)
  vhidet <- apply(vh, 3, det)
  # h <- h * ((det(H) / mean(hidet)))#^(1/d)
  vh <- sweep(vh, 3, det(H) / vhidet, "*")
  # h <- h * det(H)
  }
  
  if(opt.method == "SCALED") vh <- vh * riem.scale
  # Option 3: scale all riemannian matrix by an input constant
  
  if(missing(eval.points)) {
    # grid points for multidimensional, linear method
    gx <- numeric(0)
    for (i in 1:d) {
      gx <- c(gx, list(seq(xmin[i], xmax[i], length.out = gridsize[i])))
    }
    eval.points <- expand.grid(gx)
  
    # evaluate density on grid points
    ## TODO: CHANGE FOR NEW VOLUME DENSITY FUNCTION, CURRENTLY EVAL.POINT=DATA
    z <- NULL
    for (k in 1:n) {
      hk <- vh[,,k]
      z <-  abind::abind(z, array(mvtnorm::dmvnorm(x = eval.points, mean = x[k,], sigma = h * hk), dim = gridsize), along = d + 1)  # stack array of dimension (gridsize*gridsize) with abind
    }
    z <- rowMeans(z, dims = 2, na.rm = TRUE)
    
    f <- list(x = x, eval.points = gx, estimate = z, H = vh, r = r, gridded = TRUE)
    # if(d == 2) names(f$eval.points)[1:2] <- c("x", "y")
    
  } else {
    
    # evaluate density on given data points
    # z <- NULL
    z <- matrix(0, n, n)
    for (i in 1:n) {
      hi <- vh[,,i]
      # z <- cbind(z, mvtnorm::dmvnorm(x = eval.points, mean = x[k,], sigma = hk))  # stack vector of length n
      bindex <- which(adj_matrix[i,] <= r) # only smooth over points within the bandwidth radius ## ADD AN ARGUMENT AND CHECK DIMENSION, NEW RADIUS
      hi <- hi / mean(vh[,,bindex])
      for(j in bindex) {
        hj <- vh[,,j]
        z[i,j] <- (det(hj) / det(hi)) ^ (-0.5) * (2 * pi) ^ (-d / 2) * r ^ (-0.5) * exp( -0.5 / r * t(eval.points[i,] - eval.points[j,]) %*% solve(hi) %*% (eval.points[i,] - eval.points[j,]) ) / (pnorm(1) - pnorm(-1)) # suppK = [0,1], scale to integral to 1
      }
    }
    z <- colMeans(z, na.rm = TRUE)
    
    f <- list(x = x, eval.points = eval.points, estimate = z, H = vh, r = r, gridded = FALSE)
  }
  
  return(f)
}
