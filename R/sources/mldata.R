#' Generate random data from different mapping functions
#'
#' @param N Number of data points. By default, N = 2000.
#' @param p Dimension of meta data for mapping. Fixed at 2.
#' @param meta True manifold selected from "uniform" for uniformly distributed data \code{U(0, 1)}, 
#' "copula" for Clayton copula, or "gaussian" for Gaussian mixture from four kernels.
#' @param mapping 3-D Mapping type for the meta data. 
#'
#' @return List of six components: 3-D mapping data, 2-D meta data with true densities.
#' 
#' @export
#'
#' @examples
#' library(copula)
#' a <- mldata(N = 100, meta = "copula", mapping = "Swiss Roll")
#' library(ggplot2)
#' ggplot(as_tibble(a$data), aes(x, y)) + geom_point(aes(color = a$colors))
#' b <- mldata(meta = "gaussian", mapping = "S Curve")
#' ggplot(as_tibble(b$data), aes(x, y)) + geom_point(aes(color = b$colors))
#' dat <- as.data.frame(b$data)
#' library(scatterplot3d)
#' scatterplot3d(x = dat$x, y = dat$y, z = dat$z, color = b$colors)
mldata <- function(N = 2000, p = 2, meta = c("uniform", "copula", "gaussian"), 
                   mapping = c("Swiss Roll", "semi-sphere", "Twin Peak", "S Curve")) {
  switch (meta,
          "uniform" = {
            x = runif(N); y = runif(N)
            co <- cbind(x, y, den = 1*1)
          },
          "copula" = {
            cl2 <- claytonCopula(2, dim = p)
            co <- rCopula(N, cl2)
            co <- cbind(co, dCopula(co, cl2))
            colnames(co) <- c("x", "y", "den")
          },
          "gaussian" = {
            n <- round(N/4)
            R <- matrix(c(.02, 0,
                          0, .02), # variance-covariance matrix
                        nrow = 2, ncol = 2)
            # mu <- matrix(#~mu1, ~mu2,
            #               c(7.5, 7.5,
            #               7.5, 12.5,
            #               12.5, 7.5,
            #               12.5, 12.5), 
            #               4, 2, byrow = TRUE)
            mu <- matrix(#~mu1, ~mu2,
              c(0.25, 0.25,
                0.25, 0.75,
                0.75, 0.25,
                0.75, 0.75), 
              4, 2, byrow = TRUE)
            co <- NULL
            for(i in 1:4) {
              a <- mvtnorm::rmvnorm(n, mean = mu[i,], sigma = R)
              # True density is the mean of densities with all four cores (same R different mus)
              den <- cbind(mvtnorm::dmvnorm(a, mean = mu[1,], sigma = R),
                          mvtnorm::dmvnorm(a, mean = mu[2,], sigma = R),
                          mvtnorm::dmvnorm(a, mean = mu[3,], sigma = R),
                          mvtnorm::dmvnorm(a, mean = mu[4,], sigma = R)
              ) %>% rowMeans()
              a <- cbind(a, den, i)
              co <- rbind(co, a)
            }
            colnames(co) <- c("x", "y", "den", "label")
            co <- co %>%
              as_tibble() %>% 
              # mutate(# label = as.factor(label),
              #        x = (x - 5)/10, y = (y - 5)/10,
              #        den = den * 10 # y = h(x) = (x-5)/10; x = 10*y + 5; dx/dy = 10; f_Y(y) = f_X(x) * 10
              # ) %>% 
              as.matrix()
          }
  )
  
  switch(mapping,
         "Swiss Roll" = {
           u <- (1.5 * pi) * (1 + 2 * co[, 1:2]) 
           data <- swissRollMapping(u[,1], u[,2])
         }, 
         "semi-sphere" = {
           # u <- tibble(x = co[,1], y = co[,2]) %>% 
           #   mutate(x = 2 * pi * x, y = 2 * y - 1)
           u <- cbind(x = pi * co[,1], y = 2 * co[,2] - 1)
           data <- sphereMapping(u[,1], u[,2])
         }, 
         "Twin Peak" = {
           u <- 2 * co[, 1:2] - 1
           data <- twinPeaksMapping(u[,1], u[,2])
         }, 
         "S Curve" = {
           u <- cbind(1.5 * pi * (2 * co[,1] - 1), 2 * co[,2])
           data <- sCurveMapping(u[,1], u[,2])
         })
  
  colnames(data) <- c("x", "y", "z")
  colnames(u) <- c("x", "y") # meta data
  ifelse(meta != "gaussian", colors <- co[, "den"], colors <- as.factor(co[, "label"])) # for Gaussian kernels, use labels from four kernel indexes for coloring
  
  return(list(data = data, metadata = u, den = co[, "den"], colors = colors))
}

## Utility functions
# swiss roll mapping
swissRollMapping <- function (x, y) { # U(1.5pi, 4.5pi)
  cbind(x = x * cos(x),
        y = y,
        z = x * sin(x))
}
# semi-sphere mapping
sphereMapping <- function (phi, psi) { # y~U(-1,1), x~U(0,pi)
  cbind(x = cos(phi) * sin(psi),
        y =sin(phi) * sin(psi),
        z = cos(psi)
          )
  # cbind(x = x,   # x,y~U(-1,1)
  #       y = y, 
  #       z = sqrt(1 - (x^2 +y^2))
  # )
}
# twin peak mapping
twinPeaksMapping <- function (x, y) { # x,y~U(-1,1)
  cbind(x = x,
        y = y,
        z = sin(pi * x) * tanh(3 * y))
}
# S curve mapping
sCurveMapping <- function (t, y) { # x~U(-1.5pi, 1.5pi); y~U(0,2)
  cbind(x = sin(t),
        y = y,
        z = sign(t) * (cos(t) - 1))
}
