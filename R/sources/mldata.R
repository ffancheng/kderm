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
#' mldata(N = 1000, meta = "copula", mapping = "Swiss Roll")
#' mldata(meta = "gaussian", mapping = "S Curve")
mldata <- function(N = 2000, p = 2, meta = c("uniform", "copula", "gaussian"), 
                   mapping = c("Swiss Roll", "semi-sphere", "Twin Peak", "S Curve")) {
  switch (meta,
          "uniform" = {co = cbind(x = runif(N), y = runif(N))},
          "copula" = {
            cl2 <- claytonCopula(2, dim = p)
            co <- rCopula(N, cl2)
            co <- cbind(co, dCopula(co, cl2))
            colnames(co) <- c("x", "y", "den")
          },
          "gaussian" = {
            n <- round(N/4)
            R <- matrix(c(1, 0,
                          0, 1), 
                        nrow = 2, ncol = 2)
            # mu1 <- c(7.5, 7.5)
            # mu2 <- c(7.5, 12.5)
            # mu3 <- c(12.5, 7.5)
            # mu4 <- c(12.5, 12.5)
            mu <- tribble(
              ~mu1, ~mu2,
              7.5, 7.5,
              7.5, 12.5,
              12.5, 7.5,
              12.5, 12.5)
            co <- NULL
            for(i in 1:4) {
              # mui <- get(paste0("mu", i))
              mui <- as.matrix(mu[i,])
              a <- mvtnorm::rmvnorm(n, mean = mui, sigma = R)
              den <- mvtnorm::dmvnorm(a, mean = mui, sigma = R)
              a <- cbind(a, den, i)
              co <- rbind(co, a)
            }
            colnames(co) <- c("x", "y", "den", "label")
            co <- co %>%
              as_tibble() %>% 
              mutate(# label = as.factor(label),
                     x = (x - 5)/10, y = (y - 5)/10,
                     den = den * 10 # y = h(x) = (x-5)/10; x = 10*y + 5; dx/dy = 10; f_Y(y) = f_X(x) * 10
              ) %>% 
              as.matrix()
          }
  )
  
  switch(mapping,
         "Swiss Roll" = {
           u <- (1.5 * pi) * (1 + 2 * co) 
           data <- swissRollMapping(u[,1], u[,2])
         }, 
         "semi-sphere" = {
           u <- tibble(x = co[,1], y = co[,2]) %>% 
             mutate(x = 2 * pi * x, y = 2 * y - 1)
           data <- sphereMapping(u$x, u_$y)
         }, 
         "Twin Peak" = {
           u <- 2 * co - 1
           data <- twinPeaksMapping(u[,1], u[,2])
         }, 
         "S Curve" = {
           u <- cbind(1.5 * pi * (2 * co[,1] - 1), 2 * co[,2])
           data <- sCurveMapping(u[,1], u[,2])
         })
  
  # colors <- co[, "den"]
  ifelse(meta == "copula", colors <- co[, "den"], colors <- as.factor(co[, "label"])) # for gaussian kernels, use labels for four kernels for coloring
  
  return(list(data = data, metadata = u, colors = colors, den = co[, "den"]))
}

# swiss roll mapping
swissRollMapping <- function (x, y) { # U(1.5pi, 4.5pi)
  cbind(x = x * cos(x),
        y = y,
        z = x * sin(x))
}
# semi-sphere mapping
sphereMapping <- function (phi, psi) { # x~U(0,2pi); y~U(-1,1)
  cbind(x = cos(phi) * sin(psi),
        y = sin(phi) * sin(psi),
        z = cos(psi))
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
