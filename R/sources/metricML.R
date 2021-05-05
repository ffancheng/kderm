#' Metric manifold learning algorithm modified from the python megaman package
#'
#' @param x A set of \code{n} data points in \mathcal{R}^r
#' @param s The number of the dimensions of the embedding
#' @param k The number of nearest neighbors in manifold learning
#' @param radius The bandwidth parameter for radius nearest neighbor searching 
#' @param method The manifold learning algorithm to be applied to the data \code{x}
#' @param annmethod The approximate nearest neighbor searching method to be applied in manifold learning. The three methods are \code{"kdtree"}, \code{"annoy"}, and \code{"hnsw"}
#' @param epsilon The parameter for k-d trees to search within \code{(1+epsilon)*distance} radius
#' @param nt The number of trees parameter for Annoy algorithm
#' @param nlinks The number of links parameter for HNSW algorithm
#' @param distance The distance measure to be used in finding nearest neighbors, either \code{"euclidean"} or \code{"manhattan}
#' 
#' @return A list of the embedding coordinates \code{fn} and the embedding metric \code{hn} for each point \code{p} \in \code{x}. \code{fn} is a matrix of dimension \code{n} \times \code{n}, while \code{hn} is an array of dimension \code{n} \times \code{s} \times \code{s}
#' 
#' @references Perraul-Joncas, D., and Meila, M. (2013), "Non-linear dimensionality reduction: Riemannian metric estimation and the problem of geometric discovery," arXiv:1305.7255[stat.ML].
#' 
#' @examples
#' x <- dimRed::loadDataSet("Swiss Roll")
#' metricML(x, s = 2, k = 5, radius = .1, method = "annIsomap", annmethod = "kdtree", epsilon = 0, distance = "euclidean", treetype = "bd", searchtype = "radius")
#' metricML(x, s = 3, k = 10, method = "annLLE", annmethod = "annoy", nt = 50, distance = "manhattan")
#' 
metricML <- function(x, s, k, radius, method, annmethod = c("kdtree", "annoy", "hnsw"), epsilon = 0, nt = 50, 
                     nlinks = 16, ef.construction = 200, distance = c("euclidean", "manhattan"), diag = FALSE,
                     treetype = c("kd", "bd"), searchtype = c("standard", "priority", "radius")){
  
  # if(is.null(k) == is.null(radius)) stop("Please specify either k or radius for searching nearest neighbors, but not both. ") # ADDED IN makeKNNgraph()
  
  N <- nrow(x)
  
  ###--------------------------
  ## Step1: similarity matrix
  ###--------------------------
  nn2res <- makeKNNgraph(x = x,
                         k = k,
                         annmethod = annmethod,
                         distance = distance,
                         eps = epsilon,
                         radius = 0,
                         nt = nt,
                         nlinks = nlinks,
                         ef.construction = ef.construction,
                         treetype = treetype,
                         searchtype = searchtype,
                         )$nn2res
  
  # For radius search only, radius need to be a small value but large enough to make a connected graph
  g <- igraph::make_empty_graph(N, directed = TRUE)
  g[from = if (diag) rep(seq_len(N), times = k + 1) else rep(seq_len(N), times = k),
    to   = if (diag) as.vector(nn2res$nn.idx)  else as.vector(nn2res$nn.idx[, -1]),
    attr = "weight"] <- 
    if (diag)  as.vector(exp(- nn2res$nn.dists / radius)) else as.vector(exp(- nn2res$nn.dists[, -1] / radius)) # k_radius(p,p')
  # g <- igraph::as.undirected(g, mode = "collapse", edge.attr.comb = "first")
              
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
                     knn = k,
                     ndim = s,
                     annmethod = annmethod,
                     radius = radius, 
                     epsilon = epsilon, 
                     nt = nt, 
                     nlinks = nlinks, 
                     ef.construction = ef.construction, 
                     distance = distance,
                     treetype = treetype,
                     searchtype = searchtype,
                     .mute = c("output")
                     
  )
  fn <- e@data@data
  
  ###--------------------------
  # Step4: embedding metric hn, inverse of the Riemannian matrix, symmetric
  ###--------------------------
  hn <- riemann_metric(Y = fn, laplacian = Ln, ndim = s, invert.h = FALSE) # array of N*s*s
  
  return(list(fn_p, hn_p))
}



# Function for graph Laplacian
# Input: W: weight matrix for the neighborhood graph, radius: bandwidth parameter
# Output: L: graph Laplacian matrix
Laplacian <- function(W, radius){
  
  Tn <- Matrix::Diagonal(x=colSums(W))
  W1 <- solve(Tn) %*% W %*% solve(Tn)
  Tn1 <- Matrix::Diagonal(x=colSums(W1))
  L <- (Matrix::Diagonal(nrow(W)) - solve(Tn1) %*% W1) / radius
  
  return(L)
}



# Function for Riemannian metric for each point
riemann_metric <- function(Y, laplacian, ndim, invert.h = FALSE){
  
  H <- array(NA, dim = c(ndim, ndim, nrow(Y)))
  
  for (i in 1:ndim) {
    for (j in i:ndim) {
      yij <- Y[, i] * Y[, j]
      H[i, j, ] <- as.vector(0.5 * (laplacian %*% yij - Y[, i] * (laplacian %*% Y[, i]) - Y[, j] * (laplacian %*% Y[, j])))
      H[j, i, ] <- H[i, j, ] # symmetric matrix
    }
  }
  
  if(invert.h){
    for (i in 1:nrow(Y)) {
      H[, , i] <- MASS:ginv(H[, , i])
    }
  }
  
  return(H)
}

