#' Change the list of nearest neighbor index and distances to square distance matrix
#'
#' @param nn A list of nearest neighbor index and distances. 
#' @param sparse Whether a sparse distance matrix should be returned, TRUE by default. 
#'
#' @return
#' @export
#'
#' @examples
nn2dist <- function(nn, sparse = TRUE) {

  index <- nn$nn.idx
  N <- nrow(index)
  k <- ncol(index) - 1
  
  closest <- 
    sapply(nn, cbind) %>%
    as_tibble() %>% 
    mutate(row.idx = rep(1:N, times = k+1)) %>% 
    filter(nn.idx!=0) %>% 
    mutate(weights = nn.dists) %>%  # the weights
    arrange(row.idx)
  
  # Now construct the graph
  g <- igraph::make_empty_graph(N, directed = TRUE)
  g[from = closest$row.idx,
    to   = closest$nn.idx,
    attr = "weight"] <-
    closest$weights # k_radius(p,p')
  if(!is.connected(g)) stop("Neighborhood graph not connected. Please select a larger k/radius. ")
  
  # Distance matrix of dimension N*N
  Kn <- igraph::as_adjacency_matrix(g, attr = "weight", sparse = sparse) # dgCMatrix, or g[]
  
  return(Kn)
}
