#' Plot ellipses using the Riemannien matrix for data points in the scatterplot
#' @param x List object from `metricML()` function containing embeddings and Riemannian metric
#' @param add Whether the ellipses are to be added to another ggplot without the scatterplot. 
#' By default `add=FALSE` generates a scatterplot with some ellipses.
#' @param n.plot Number of randomly selected data points to add ellipses
#' @param ell_size Scaling parameter for the size of the ellipses
#' @param color Color of the ellipses
#' @param fill Fill color of the ellipses
#' @param alpha Transparency for the ellipses
#' @param ... 
#'
#' @return A ggplot object of a 2-d scatterplot with ellipses added to some points
#' @export
#'
#' @examples
#' library(ggplot2)
#' library(dimRed)
#' x <- dimRed::loadDataSet("Swiss Roll")
#' a <- metricML(x, s = 2, k = 5, radius = 10, method = "annIsomap", annmethod = "kdtree", epsilon = 0, distance = "euclidean", treetype = "kd", searchtype = "radius")
#' plot_ellipse(a)
#' plot_embedding(a) + plot_ellipse(a, add = TRUE)
plot_ellipse <- function(x, add = FALSE, ell.no = 50, ell.size = 1,
             color = blues9[5], fill = blues9[5], alpha = 0.2, ...){
  fn <- x$embedding
  colnames(fn) <- paste0("E", 1:ncol(fn))
  stopifnot(ncol(fn)==2)
  Rn <- x$rmetric # array
  e <- riem2ellipse(Rn, ell.size) %>% 
    cbind(fn) %>% 
    as_tibble()
  
  p <- geom_ellipse(data = slice_sample(e, n = ell.no), 
                    aes(x0 = E1, y0 = E2, a = a, b = b, angle = angle, group = -1), 
                    color = color, fill = fill, alpha = alpha, ...)
  if(!add) {
    p <- plot_embedding(fn, embedding = TRUE) + p
  }
 
  return(p)
}
