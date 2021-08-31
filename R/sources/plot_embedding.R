# Embedding plot of metricML() output for electricity data
plot_embedding <- function(x, embedding = FALSE, color = NULL) {
  
  if(embedding){
    fn <- x
  } else{
    fn <- x$embedding
  }
  N <- nrow(fn)
  
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
    ggplot(aes(x = E1, y = E2, col = color)) + 
    geom_point() + 
    coord_fixed(ratio = 1) + 
    # color + 
    theme(legend.position = 'bottom')
  
  return(p)
}
