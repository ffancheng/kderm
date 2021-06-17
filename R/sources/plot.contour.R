# takes the output of metricML() as input for contour plot
# n.grid for grid size
plot.contour <- function(x, n.grid = 50){
  
  fn <- x$embedding
  Rn <- x$rmetric # array
  # h <- t(apply(Rn, 3, diag))
  f <- vkde2d(x = fn[,1], y = fn[,2], h = Rn, n = n.grid)
  # str(f)
  # image(f)
  p <- filled.contour(f, color.palette = viridis,
                      plot.axes = { axis(1); axis(2);
                        points(fn, pch = 3, col= hcl(c=20, l = 8))})
}



# use akima for interpolation
# plot.contour <- function(x, f, pixel=100, lenbreak=5, plot.hdr = FALSE, ...){
#   embed_den <- as_tibble(cbind(x = x[,1], y = x[,2], z = f)) %>%
#     drop_na()
#   pixel <- 100
#   lenbreak <- 5
#   akima.spl <- akima::interp(embed_den$x, embed_den$y, embed_den$z, nx=pixel, ny=pixel, linear=FALSE) # akima uses splines for interpolation, which is only useful for area with enough data points, but is not for points on the edges/no points
#   
#   p1 <- NULL
#   if(plot.hdr){
#     p1 <- hdrscatterplot(embed_den$x, embed_den$y, noutliers = 10)
#   }
#   
#   p <- filled.contour(akima.spl, color.palette = viridis,
#                  plot.axes = { axis(1); axis(2);
#                    title("smooth interp(*, linear = FALSE)");
#                    points(embed_den, pch = 3, col= hcl(c=20, l = 10))}, 
#                  ...)
#   
#   return(p)
# }

# x <- dimRed::loadDataSet("Swiss Roll")
# N <- nrow(x)
# s <- 2
# k <- 20
# method <- "annIsomap"
# annmethod <- "kdtree"
# distance <- "euclidean"
# treetype <- "kd"
# searchtype <- "radius" # change searchtype for radius search based on `radius`, or KNN search based on `k`
# radius <- .4 # the bandwidth parameter, \sqrt(\elsilon), as in algorithm
# metric_isomap <- metricML(x, s = s, k = k, radius = radius, method = method, annmethod = annmethod, eps = 0, distance = distance, treetype = treetype, searchtype = searchtype, invert.h = TRUE)
# f <- mkde(x = metric_isomap$embedding, h = metric_isomap$rmetric)
# plot.contour(x = metric_isomap$embedding, f, plot.hdr = F)

