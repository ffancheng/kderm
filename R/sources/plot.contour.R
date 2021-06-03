plot.contour <- function(x, f, pixel=100, lenbreak=5, plot.hdr = FALSE, ...){
  embed_den <- as_tibble(cbind(x = x[,1], y = x[,2], z = f)) %>%
    drop_na()
  pixel <- 100
  lenbreak <- 5
  akima.spl <- interp(embed_den$x, embed_den$y, embed_den$z, nx=pixel, ny=pixel, linear=FALSE)
  
  p1 <- NULL
  if(plot.hdr){
    p1 <- hdrscatterplot(embed_den$x, embed_den$y, noutliers = 10)
  }
  
  p <- filled.contour(akima.spl, color.palette = viridis,
                 plot.axes = { axis(1); axis(2);
                   title("smooth interp(*, linear = FALSE)");
                   points(embed_den, pch = 3, col= hcl(c=20, l = 10))}, 
                 ...)
  
  return(list(p=p))
}
