# Variable kernel density estimate
# \hat{f}(x) = n^{-1} \sum_i K_H (x - X_i)

# calculate the density estimate for each point in emb_isomap with each H matrix, H_isomap[i,,]

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
      s[i] <- s[i] + abs(det(h[,,i]))^(-1/2) * exp(- 1/2 * t(x[i,] - x[j,]) %*% solve(h[,,i]) %*% (x[i,] - x[j,]))
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
