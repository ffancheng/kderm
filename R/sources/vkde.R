# input x is a multivariate data/matrix
# takes care of the computation when bandwidth is full matrix instead of fixed vector or diagonal matrix
# lims is a 2-column matrix of ranges(xmin, xmax) for each column in x
vkde <- function(x, h, gridsize = 25, xmin = apply(x, 2, min), xmax = apply(x, 2, max), eval.points = NULL){

  nx <- nrow(x)
  d <- ncol(x)
  if(any(!is.finite(x)))
    stop("missing or infinite values in the data are not allowed")
  if(any(!is.finite(lims)))
    stop("only finite values are allowed in the minimum/maximum values for grids")

  if(is.vector(h)) return(ks::kde(x, h = h, gridsize = gridsize, xmin = xmin, xmax = xmax, eval.points = eval.points)) # fixed diagonal bandwidth # return(MASS::kde2d(x, y, h, n, lims)) # only works for 2d

  # TODO: grid points for multidimension
  nd <- rep(gridsize, length.out = d)  # number of grids for all dimensions
  gx <- seq.int(lims[1L], lims[2L], length.out = nd[1L])
  gy <- seq.int(lims[3L], lims[4L], length.out = nd[2L])
  g <- expand.grid(gx, gy)

  a <- cbind(xmin, xmax) %>%
    as_tibble() %>%
    mutate(list(seq.int(xmin, xmax, length.out = d)))

  g <- seq.int(xmin[1], xmax[1], length.out = nd[1])
  for(i in 2:d){
    gy <- seq.int(xmin[i], xmax[i], length.out = nd[i])
    g <- expand.grid(g, gy)
  }

  g <- seq.int(xmin, xmax, length.out = gridsize)



  # bandwidth matrix is a 2*1?
  # h <- if (missing(h)) t(replicate(nx, c(bandwidth.nrd(x), bandwidth.nrd(y)))) else h
  # else rep(h, length.out = 2L) # for fixed h input in kde2d

  # if (any(h <= 0))
  #   stop("bandwidths must be strictly positive")

  # # if bandwidth is a d*d*nx array
  z <- array(NA, c(n[1L], n[2L], nx))

  for (k in 1:nx) {
    hk <- h[,,k]
    z[,,k] <- det(hk) ^ (-1/2) * matrix(mvtnorm::dmvnorm(x = g, mean = x[k,], sigma = hk), nrow = nd[1], byrow = FALSE) # scalar
  }

  z <- rowMeans(z, dims = 2, na.rm = TRUE)

  ## TODO: optimize hk as bandwidth pointwise
  # AMISE



hkd <- 1 / (2 * pi) * det(h[,,k])^(-1/2)
hks <- solve(h[,,k])
xyk <- c(x[k], y[k])

for(i in 1:n[1L]) {
  for (j in 1:n[2L]) {

    gxy <- c(gx[i], gy[j])
    s[i, j, k] <- hkd * exp( - 1/2 * t(gxy - xyk) %*% hks %*% (gxy - xyk))
  }
}
}




# z <- rowMeans(s, dims = 2)

# for(i in 1:n[1L]) {
#   for (j in 1:n[2L]) {
#
#     s = 0
#     for (k in 1:nx) {
#       # s <- s + 1 / (2 * pi) * (-1/(2*h[k,1]*h[k,2])) * exp(- 1/2 * t(gx[i] - x[k]) * 1/(h[k, 1]*h[k, 2]) * (gy[j] - y[k]))
#       s <- s + 1 / (2 * pi) * det(h[,,k])^(-1/2) *
#         exp( - 1/2 * t(c(gx[i], gy[j]) - c(x[k], y[k])) %*% solve(h[,,k]) %*% ( c(gx[i], gy[j]) - c(x[k], y[k]) ) )
#     }
#     z[i, j] <- s / nx
#
#   }
# }

# from MASS:: kde2d()
# h <- h/4                            # for S's bandwidth scale
# ax <- outer(gx, x, "-" )/h[1L]
# ay <- outer(gy, y, "-" )/h[2L]
# z <- tcrossprod( matrix(dnorm(ax), , nx)/h[1L], matrix(dnorm(ay), , nx)/h[2L] ) /
#   nx
# # z <- tcrossprod(matrix(dnorm(ax), , nx), matrix(dnorm(ay), , nx)) /
# #   (nx * h[1L] * h[2L])

#   return(list(x = gx, y = gy, z = z, grid = g))
# }

