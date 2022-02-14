# Convert a 2*2 covariance matrix to an ellipse's radii lengths and angle
cov2ellipse <- function(x) {
  aa <- x[1, 1]
  bb <- x[1, 2]
  cc <- x[2, 2]
  a2 <- (aa + cc) / 2 + sqrt(((aa - cc) / 2) ^ 2 + bb ^ 2) 
  b2 <- (aa + cc) / 2 - sqrt(((aa - cc) / 2) ^ 2 + bb ^ 2) 
  angle <- atan((a2 - aa) / bb)  # Or: atan(bb / (a2 - aa))
  # radius <- sqrt(qchisq(p=0.95, df = 2)) # radius to cover 95% of all data, error/confidence ellipse
  # radius <- 1
  
  return(list(a = sqrt(a2),# * radius, 
              b = sqrt(b2),# * radius, 
              angle = angle))
}

# Take Riemmanian metric from an array to a tibble for ggforce::geom_ellipse()
riem2ellipse <- function(x, ell.size = 1){
  Rn <- purrr::array_branch(x, 3)
  # convert list of covariance matrix to tibble of a, b, angle
  e <- sapply(Rn, cov2ellipse) %>% 
    t() %>%
    apply(2, unlist) %>% 
    # cbind(tod = rep(1:48, N / 48)) %>% 
    as_tibble() %>% 
    mutate(a = a * ell.size,
           b = b * ell.size)
  
  return(e)
}

# library(tibble)
# library(MASS)
# library(ggforce)
# sigma=matrix(c(1.4,-0.5,-.5,0.4), nrow = 2, ncol = 2)
# samples<-mvrnorm(n=100, c(1,2), sigma)
# colnames(samples)<-c("X","Y")
# samples_df = as_tibble(samples)
# 
# my_cov<-cov(samples_df)
# par_elipsy<-cov2elipse(my_cov)
# 
# ggplot(samples_df, aes(x=X, y=Y)) +
#   geom_point() +
#   geom_ellipse(mapping = aes(x0=1, y0=2, a=par_elipsy$a, b=par_elipsy$b, angle=par_elipsy$angle), color="blue", fill = blues9[3], alpha = 0.01)
