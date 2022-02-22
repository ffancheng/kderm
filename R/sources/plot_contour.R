# takes the output of metricML() as input for contour plot
# n.grid for grid size
#' Title
#'
#' @param x List object from `metricML()` function containing the 2-d embeddings and the Riemannian metric
#' @param gridsize Grid size for estimating densities
#' @param f Pre-computed densities
#' @param riem.scale Scaling parameter for all Riemannien metric in variable KDE
#'
#' @return A full page display of the filled contour with scatterplot of embeddings
#' @export
#'
#' @examples
#' 
plot_contour <- function(x, gridsize = 20, f = NULL, riem.scale = 1){
  
  fn <- x$embedding
  Rn <- x$rmetric # array
  # h <- t(apply(Rn, 3, diag))
  
  if(is.null(f)) f <- vkde(x = fn, h = Rn * riem.scale, gridsize = gridsize)
  # f <- vkde2d(x = fn[,1], y = fn[,2], h = Rn * riem.scale, gridsize = gridsize)
  
  xyz <- list(x = f$eval.points[[1]], y = f$eval.points[[2]], z = f$estimate)
  # image(f)
  filled.contour(xyz, color.palette = viridis,
                 plot.axes = { axis(1); axis(2); points(fn, pch = 3, col= hcl(c=20, l = 8, alpha=0.2))}
  )
  
}
