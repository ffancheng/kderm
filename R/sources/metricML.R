#' Metric manifold learning algorithm modified from the python megaman package
#'
#' @param x A set of \code{n} data points in \mathcal{R}^r
#' @param s The number of the dimensions of the embedding
#' @param k The number of nearest neighbors in manifold learning
#' @param radius The bandwidth parameter for radius nearest neighbor searching 
#' @param method The manifold learning algorithm to be applied to the data \code{x}
#' @param annmethod The approximate nearest neighbor searching method to be applied in manifold learning. The three methods are \code{"kdtree"}, \code{"annoy"}, and \code{"hnsw"}, but only \code{"kdtree"} is implemented for radius search now
#' @param eps The parameter for k-d trees to search within \code{(1+epsilon)*distance} radius
#' @param nt The number of trees parameter for Annoy algorithm
#' @param nlinks The number of links parameter for HNSW algorithm
#' @param distance The distance measure to be used in finding nearest neighbors, either \code{"euclidean"} or \code{"manhattan}
#' @param treetype Character vector specifying the standard 'kd' tree or a 'bd' (box-decomposition, AMNSW98) tree which may perform better for larger point sets. See details in ?RANN::nn2
#' @param searchtype Search types: priority visits cells in increasing order of distance from the query point, and hence, should converge more rapidly on the true nearest neighbour, but standard is usually faster for exact searches. radius only searches for neighbours within a specified radius of the point. If there are no neighbours then nn.idx will contain 0 and nn.dists will contain 1.340781e+154 for that point. See details in ?RANN::nn2
#' 
#' @return A list of the embedding coordinates \code{fn} and the embedding metric \code{hn} for each point \code{p} \in \code{x}. \code{fn} is a matrix of dimension \code{n} \times \code{n}, while \code{hn} is an array of dimension \code{n} \times \code{s} \times \code{s}
#' 
#' @references Perraul-Joncas, D., and Meila, M. (2013), "Non-linear dimensionality reduction: Riemannian metric estimation and the problem of geometric discovery," arXiv:1305.7255[stat.ML].
#' 
#' @examples
#' require(dimRed)
#' x <- dimRed::loadDataSet("Swiss Roll")
#' metricML(x, s = 2, k = 5, radius = .4, method = "annIsomap", annmethod = "kdtree", epsilon = 0, distance = "euclidean", treetype = "kd", searchtype = "radius")
#' metricML(x, s = 3, k = 10, method = "annLLE", annmethod = "annoy", nt = 50, distance = "manhattan")
#' 
metricML <- function(x, s, k = min(10, nrow(x)), radius = 0, 
                     # adjacency_graph = NULL, 
                     method, annmethod = c("kdtree", "annoy", "hnsw"), 
                     eps = 0, nt = 50, nlinks = 16, ef.construction = 200,
                     distance = c("euclidean", "manhattan"), diag = FALSE,
                     treetype = c("kd", "bd"),
                     searchtype = c("standard", "priority", "radius"),
                     ...
                     ){
  
  N <- nrow(x)
  
  # if(is.null(k) == is.null(radius)) stop("Please specify either k or radius for k-d trees to find nearest neighbors, but not both. ")
  # When there are more than `k` NNs found using radius search, print k NNs instead of 10. RANN::nn2() only prints `k` NNs
  if(searchtype == "radius"){
    if(radius <= 0) stop("Please specify a positive value for `radius` when using radius search. ")
    k <- N - 1
  }
  
  # TODO: input as the NN graph, skip Step1, but assign weights
  # if(is.null(adjacency_graph)){
  #
  # }
  
  ###--------------------------
  ## Step1: similarity matrix, symmetric
  ###--------------------------
  nn2res <- dplyr::case_when(distance=="euclidean" ~ 
                               RANN::nn2(data = x, query = x, k = k + 1, 
                                         treetype = treetype, searchtype = searchtype, eps = eps,
                                         radius = radius),
                             distance=="manhattan" ~ 
                               RANN.L1::nn2(data = x, query = x, k = k + 1, 
                                            treetype = treetype, searchtype = searchtype, eps = eps,
                                            radius = radius),
  )
  names(nn2res) <- c("nn.idx", "nn.dists")
  
  # Convert RANN::nn2() output as N*N adjacency matrix
  # Kn <- Matrix::Matrix(0, N, N, sparse = TRUE)
  closest <- 
    sapply(nn2res, cbind) %>%
    as_tibble() %>% 
    mutate(row.idx = rep(1:N, times = k+1)) %>% 
    filter(nn.idx!=0) %>% 
    mutate(weights = exp(- nn.dists / (radius^2))) %>% 
    arrange(row.idx)
  # View(closest)
  
  # Now construct the graph
  g <- igraph::make_empty_graph(N, directed = TRUE)
  g[from = closest$row.idx,
    to   = closest$nn.idx,
    attr = "weight"] <-
    closest$weights # k_radius(p,p')
  # is.connected(g)
  
  # for(i in 1:nrow(closest)){
  #   Kn[closest$row.idx[i], closest$nn.idx[i]] <- closest$weights[i]
  # }
  # # Kn[1:10, 1:10]
  # # isSymmetric(Kn) # TRUE
  
  # # TODO: The following code is for creating neighborhood graph, but errors apprear for radius search when there are 0s in the nn.idx from RANN::nn2()
  # # For radius search only, radius need to be a small value but large enough to make a connected graph
  # g <- igraph::make_empty_graph(N, directed = TRUE)
  # g[from = if (diag) rep(seq_len(N), times = k + 1) else rep(seq_len(N), times = k),
  #   to   = if (diag) as.vector(nn2res$nn.idx)  else as.vector(nn2res$nn.idx[, -1]),
  #   attr = "weight"] <- 
  #   if (diag)  as.vector(exp(- nn2res$nn.dists / (radius^2))) else as.vector(exp(- nn2res$nn.dists[, -1] / (radius^2))) # k_radius(p,p')
  # 
  # # g <- igraph::as.undirected(g, mode = "collapse", edge.attr.comb = "first")
  #             
  Kn <- igraph::as_adjacency_matrix(g, attr = "weight", sparse = TRUE) # dgCMatrix, or g[]
  
  
  ###--------------------------
  ## Step2: Laplacian matrix
  ###--------------------------
  # igraph::laplacian_matrix(g)
  Ln <- Laplacian(W = Kn, radius = radius)
  

  ###--------------------------
  ## Step3: embedding coordinates fn
  ###--------------------------
  e <- dimRed::embed(x,
                     .method = method,
                     knn = 100,
                     ndim = s,
                     annmethod = annmethod,
                     radius = radius, 
                     eps = eps,
                     nt = nt,
                     nlinks = nlinks,
                     ef.construction = ef.construction,
                     distance = distance,
                     treetype = treetype,
                     searchtype = searchtype,
                     .mute = c("output"),
                     ...
                     
  )
  fn <- e@data@data
  colnames(fn) <- paste0("E", 1:s)
  
  ###--------------------------
  # Step4: embedding metric hn, inverse of the Riemannian matrix, symmetric
  ###--------------------------
  hn <- riemann_metric(Y = fn, laplacian = Ln, ndim = s, invert.h = T) # array of N*s*s
  
  return(list(embedding=fn, rmetric=hn, 
              weighted_graph=g,
              adj_matrix=Kn))
}



# Function for graph Laplacian
# Input: W: N*N weight matrix for the neighborhood graph, radius: bandwidth parameter
# Output: L: N*N graph Laplacian matrix
Laplacian <- function(W, radius){
  
  Tn <- Matrix::Diagonal(x = 1 / rowSums(W)) # inverse of a diagonal matrix
  # temp <- solve(Tn)
  W1 <- Tn %*% W %*% Tn
  Tn1 <- Matrix::Diagonal(x = 1 / rowSums(W1)) # inverse of Tn1
  L <- (Matrix::Diagonal(nrow(W)) - Tn1 %*% W1) / radius
  
  return(L)
}



# Function for Riemannian metric for each point
# The Riemannian metric and its dual associated with an embedding Y. 
# The Riemannian metric is currently denoted by G, its dual by H, and the Laplacian by L. 
# G at each point is the matrix inverse of H.
riemann_metric <- function(Y, laplacian, ndim, invert.h = TRUE){
  
  H <- array(NA, dim = c(ndim, ndim, nrow(Y)))
  
  for (i in 1:ndim) {
    for (j in i:ndim) {
      yij <- Y[, i] * Y[, j]
      H[i, j, ] <- as.vector(0.5 * (laplacian %*% yij - Y[, i] * (laplacian %*% Y[, i]) - Y[, j] * (laplacian %*% Y[, j])))
      H[j, i, ] <- H[i, j, ] # symmetric matrix
    }
  }
  
  # Array H corresponds to \tilde{H} in Step 4(a)
  # The embedding metric H is the pseudo inverse of \tilde{H}
  # for (i in 1:nrow(Y)) {
  #   Hsvals <- eigen(H[,,i])
  #   Huu <- Hsvals$vectors
  #   Hvv <- Hsvals$values[1:ndim]
  #   Hvv1 <- diag(x = 1 / Hvv)
  #   H[,,i] <- Huu %*% Hvv1 %*% t(Huu)
  #   H[, , i] <- 0.5 * (H[, , i] + t(H[, , i]))
  # }
  
  
  if(invert.h){
    for (i in 1:nrow(Y)) {
      H[, , i] <- solve(H[, , i])
      H[, , i] <- 0.5 * (H[, , i] + t(H[, , i])) # fix not symmetric issue
    }
  }

  return(H)
}

# a <- H[,,1]
# eigen(a)
# is.positive.definite(a) # FALSE
