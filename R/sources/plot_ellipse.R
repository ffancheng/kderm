plot.ellipse <- function(x, add = FALSE, n.plot = 50, scale = 20,
             color = blues9[5], fill = blues9[5], alpha = 0.2, ...){
  
  fn <- x$embedding
  Rn <- x$rmetric # array
  e <- riem2ellipse(Rn, scale = 20) %>% 
    cbind(fn) %>% 
    as_tibble()
  
  p <- geom_ellipse(data = slice_sample(e, n = n.plot), 
                    aes(x0 = E1, y0 = E2, a = a, b = b, angle = angle, group = -1), 
                    color = color, fill = fill, alpha = alpha, ...)
  
  if(add){
    p <- p
  } else {
    p <- plot.embedding(x, embedding = F) + p
  }
 
  return(p)
}
