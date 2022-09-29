# Categorize density to `length(prob)` levels prob and plot the embedding colored by the density levels
plot_embedding_color_levels <- function(x = trueden, N, prob = c(1, 50, 95, 99), noutliers, fn = NULL) {
  if(max(prob) > 50){ # Assume prob is coverage percentage
    alpha <- (100 - prob) / 100
  } else {# prob is tail probability (for backward compatibility)
    alpha <- prob
  }
  alpha <- sort(alpha) # 0.99 0.50 0.05 0.01
  falpha <- quantile(x, alpha) # increasing
  region <- numeric(N) + 100 # rep(100, N)
  for (i in seq_along(prob)) region[x > falpha[i]] <- 100 - prob[i]
  if (is.null(noutliers)) noutliers <- sum(region > max(prob))
  noutliers <- min(noutliers, N)
  
  # Construct region factor
  prob <- sort(unique(region[region < 100]), decreasing = TRUE)
  prob <- c(prob, 100)
  Region <- factor(region,
                   levels = prob,
                   labels = c(paste(head(prob, -1)), ">99") # "99"  "50"  "5"   "1"   ">99"
  )
  # Sort data so the larger regions go first (other than outliers)
  # rr <- region
  # rr[region == 100] <- 0
  # ord <- order(rr, decreasing = TRUE)
  ord <- order(x, decreasing = T)  # outliers at the end
  # fn_dec <- fn[ord,] # data is reordered by ord
  outliers <- tail(ord, noutliers)
  
  # if (noutliers > 0) {
  #   outlier_rank <- order(x[ord])  # order density as well to match the new data
  #   outliers <- outlier_rank[seq_len(noutliers)] # take top noutliers labels for annotation
  # }
  
  p <- fn %>% 
    as_tibble() %>% 
    ggplot(aes(x = E1, y = E2)) + 
    ggplot2::geom_point(ggplot2::aes_string(col = "Region")) + 
    ggplot2::scale_colour_manual(
      name = "HDRs",
      breaks = c(paste(head(sort(prob), -1)), ">99"),
      values = c(RColorBrewer::brewer.pal(length(prob), "YlOrRd")[-1], "#000000")) + 
    ggplot2::annotate("text",
                      x = fn[outliers, 1] + diff(range(fn[outliers, 1])) / 50, 
                      y = fn[outliers, 2] + diff(range(fn[outliers, 2])) / 50,
                      label = outliers, # seq_len(N)[outliers], 
                      col = "blue", cex = 2.5
    )
  return(p)
}
