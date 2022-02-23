# f: estimates for all data points x, not grid points
plot_outlier <- function(x, gridsize = 20, f = NULL, prob = c(1, 50, 99), noutliers = 20, label = NULL, riem.scale = 1, ell.size = 1, ...){
  
  fn <- x$embedding
  Rn <- x$rmetric
  # if(is.null(f)) f <- vkde2d(x = fn[,1], y = fn[,2], h = Rn * riem.scale, gridsize = gridsize) # works for 2d, use linear interpolation
  if(is.null(f)) f <- vkde(x = fn, h = Rn*riem.scale, gridsize = gridsize, eval.points = fn)
  
  den <- f
  E1 <- fn[,1]
  E2 <- fn[,2]
  # Convert prob to coverage percentage if necessary
  if(max(prob) > 50) {# Assume prob is coverage percentage
    alpha <- (100 - prob) / 100
  } else {# prob is tail probability (for backward compatibility)
    alpha <- prob}
  alpha <- sort(alpha)
  # Calculates falpha needed to compute HDR of bivariate density den.
  
  # # Also finds approximate mode.
  # fxy <- hdrcde:::interp.2d(den$x,den$y,den$z,E1,E2) # where to estimate density, stored in den
  fxy <- f$estimate # change eval.points = data in vkde2d
  falpha <- quantile(fxy, alpha)
  index <- which.max(fxy)
  
  mode <- c(E1[index],E2[index]) # maximam density point
  hdr2d_info <- structure(list(mode=mode, falpha=falpha, fxy=fxy, den=den, alpha=alpha, x=E1, y=E2), class="hdr2d") # list for hdr.2d() output
  # plot.hdr2d(hdr2d_info, show.points = T, outside.points = T, pointcol = grey(0.5), xlim = round(range(E1)), ylim = round(range(E2)))
  
  p_hdr <- hdrscatterplot_new(E1, E2, levels = prob, noutliers = noutliers, label = label, den = hdr2d_info)
  p_outlier_vkde <- p_hdr$p + 
    plot_ellipse(x, add = T, ell.size = ell.size, ...)
                 # color = blues9[5], fill = blues9[5], alpha = 0.2, ...)
  
  return(list(p = p_outlier_vkde, outlier = p_hdr$outlier, densities = fxy, hdr2d_info = hdr2d_info))
}

# Example
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
