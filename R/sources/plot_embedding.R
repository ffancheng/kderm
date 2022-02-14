### could be further implemented as plot for S3 class
#' 2-d embedding plot of metricML() output
#'
#' @param x List object from `metricML()` function containing the 2-d embeddings
#' @param embedding Whether the input is the 2-d embedding data
#' @param color Colors for each data point 
#' @param alpha Transparency for each data point
#'
#' @return A ggplot object of the scatterplot
#' @export
#'
#' @examples
#' library(dimRed)
#' x <- dimRed::loadDataSet("Swiss Roll")
#' a <- metricML(x, s = 2, k = 5, radius = .4, method = "annIsomap", annmethod = "kdtree", epsilon = 0, distance = "euclidean", treetype = "kd", searchtype = "radius")
#' plot_embedding(a)
#' plot_embedding(a$embedding, embedding = TRUE)
plot_embedding <- function(x, embedding = FALSE, color = NULL, alpha = NULL) {
  
  if(embedding){
    fn <- x
  } else{
    fn <- x$embedding
  }
  stopifnot(ncol(fn)==2)
  N <- nrow(fn)
  if(!is.null(color) & length(color)!=N)
    stop("Given colors should have the same length as the number of points!")
  
  # color for electricity data
  # if(color == NULL) {
  #   color <- colorspace::scale_color_continuous_sequential(
  #     palette = "viridis",
  #     breaks = c(12, 24, 36, 48),
  #     labels=c("06:00", "12:00", "18:00", "24:00"),
  #     name="Time of day",
  #     guide=guide_colorbar(barwidth = 10))
  # } 
  
  p <- fn %>% 
    # cbind(tod = rep(1:48, times = N / 48)) %>% 
    as_tibble() %>% 
    ggplot(aes(x = E1, y = E2, color = color, alpha = alpha)) + 
    geom_point() + 
    coord_fixed(ratio = 1) + 
    # color + 
    theme(legend.position = 'bottom')
  
  return(p)
}
