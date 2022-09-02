#' Learn Metric algorithm modified from the python megaman package
#'
#' @param x A set of \code{n} data points in \mathcal{R}^r
#' @param s The number of the dimensions of the embedding
#' @param k The number of nearest neighbors in manifold learning
#' @param radius The radius for nearest neighbor searching 
#' @param bandwidth The bandwidth parameter of the kernel for weighted graph matrix, \code{\sqrt{0.4} for heat kernel}, could be the same as \code{radius}
#' @param const The constant term for the estimeted Laplacian matrix that depends on the choice of kernel, const = 0.25 as suggested in the Learn Metric algorithm
#' @param adjacency The adjacency matrix of the data of dimension \code{n*n}, NULL by default
#' @param affinity The weighted adjacency matrix with gaussian kernel of dimension \code{n*n}, NULL by default
#' @param method The manifold learning algorithm to be applied to the data \code{x}
#' @param fn The low-dimensional embedding from manifold learning, NULL by default
#' @param annmethod The approximate nearest neighbor searching method to be applied in manifold learning. The three methods are \code{"kdtree"}, \code{"annoy"}, and \code{"hnsw"}, but only \code{"kdtree"} is implemented for radius search now
#' @param eps The parameter for k-d trees to search within \code{(1+epsilon)*distance} radius
#' @param nt The number of trees parameter for Annoy algorithm
#' @param nlinks The number of links parameter for HNSW algorithm
#' @param distance The distance measure to be used in finding nearest neighbors, either \code{"euclidean"} or \code{"manhattan}
#' @param treetype Character vector specifying the standard 'kd' tree or a 'bd' (box-decomposition, AMNSW98) tree which may perform better for larger point sets. See details in ?RANN::nn2
#' @param searchtype Search types: priority visits cells in increasing order of distance from the query point, and hence, should converge more rapidly on the true nearest neighbour, but standard is usually faster for exact searches. Radius only searches for neighbours within a specified radius of the point. If there are no neighbours then nn.idx will contain 0 and nn.dists will contain 1.340781e+154 for that point. See details in ?RANN::nn2
#' @param invert.h Whether the Riemannian metric needs to be inverted. By default it is FALSE
#' 
#' @return A list of the embedding coordinates \code{fn} and the embedding metric \code{hn} for each point \code{p} \in \code{x}. \code{fn} is a matrix of dimension \code{n} \times \code{s}, while \code{hn} is an array of dimension \code{n} \times \code{s} \times \code{s}
#' 
#' @references Perraul-Joncas, D., and Meila, M. (2013), "Non-linear dimensionality reduction: Riemannian metric estimation and the problem of geometric discovery," arXiv:1305.7255[stat.ML].
#' 
#' @examples
#' require(dimRed)
#' x <- dimRed::loadDataSet("Swiss Roll")
#' metricML(x, s = 2, k = 5, radius = .4, method = "annIsomap", annmethod = "kdtree", epsilon = 0, distance = "euclidean", treetype = "kd", searchtype = "radius")
#' metricML(x, s = 3, k = 10, method = "annLLE", annmethod = "annoy", nt = 50, distance = "manhattan")
#' 
metricML <- function(x, s = 2, k = min(10, nrow(x)), radius = 0, 
                     bandwidth = 0.4, const = 0.25, 
                     adjacency = NULL, affinity = NULL,
                     method, fn = NULL,
                     annmethod = c("kdtree", "annoy", "hnsw"),
                     eps = 0, nt = 50, nlinks = 16, ef.construction = 200, ef.search = 10,
                     distance = c("euclidean", "manhattan")[1], diag = FALSE,
                     treetype = c("kd", "bd")[1],
                     searchtype = c("standard", "priority", "radius")[1],
                     perplexity = round(k/3), theta = 0.5, # t-SNE
                     invert.h = TRUE,
                     ...
                     ) {
  
  # input as the adjacency/affinity matrix, skip Step1
  if(is.null(x)){
    
    if(!is.null(affinity)){
      if(nrow(affinity) != ncol(affinity)) stop("The affinity matrix of the data should be a square matrix. ")
      Kn <- affinity
    } 
    
    if(!is.null(adjacency)){
      if(nrow(adjacency) != ncol(adjacency)) stop("The adjacency matrix of the data should be a square matrix. ")
      affinity <- exp(- (adjacency / radius) ^ 2)
      Kn <- affinity
    } 
    
    x <- Kn
    if(is.null(affinity) & is.null(adjacency)) stop("Please specify one of the data matrix, the adjacency matrix, and the affinity matrix along with the searching radius as the inputs.")
    
  } else { # use data matrix input

    # if(is.null(k) == is.null(radius)) stop("Please specify either k or radius for k-d trees to find nearest neighbors, but not both. ")
    # When there are more than `k` NNs found using radius search, print k NNs instead of 10. RANN::nn2() only prints `k` NNs
    # if(searchtype == "radius") {
    #   if(radius <= 0) stop("Please specify a positive value for `radius` when using radius search. ")
    # }
      
    ##--------------------------
    # Step3: embedding coordinates fn
    ##--------------------------
    if(is.null(fn)) { # skipped if the embedding fn is given
      e <- embed(x,
                 .method = method,
                 knn = k,
                 ndim = s,
                 annmethod = annmethod,
                 radius = radius,
                 eps = eps,
                 nt = nt,
                 nlinks = nlinks,
                 ef.construction = ef.construction,
                 ef.search = ef.search,
                 distance = distance,
                 treetype = treetype,
                 searchtype = searchtype,
                 perplexity = perplexity,
                 theta = theta,
                 .mute = c("output"),
                 # ...
                 
      )
      fn <- e@data@data
    }
    
    # geodist <- igraph::distances(g, algorithm = "dijkstra")
    # fn <- cmdscale(geodist, k = 2)
    s <- ncol(fn)
    colnames(fn) <- paste0("E", 1:s) 
    
    N <- nrow(x)
    if (searchtype == "radius") {  
      k <- N - 1 # for printing full distance matrix, but will cause error in makeKNNgraph() when building weighted graph edges; in nn2(), k is the maximum number of NNs to output, so k is set as N even if radius is large
    }
    
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
    
    W <- matrix(0, N, N)
    for(i in 1:N) {
      W[i, nn2res$nn.idx[i,]] <- exp(-nn2res$nn.dists[i, nn2res$nn.idx[i,] != 0] / (bandwidth ^ 2))
    }

    # get pairwise distance matrix
    Kn <- nn2dist(nn2res, N = N) # TODO: check the weight matrix
    # W <- exp(- Kn / (bandwidth ^ 2)) # heat kernel for weighted graph
    Kn[Kn == 1e+05] <- 0
    # g <- igraph::graph_from_adjacency_matrix(Kn)
    
  }
      
    
  ###--------------------------
  ## Step2: Laplacian matrix
  ###--------------------------
  # igraph::laplacian_matrix(g)
  Ln <- Laplacian(W = W, bandwidth = bandwidth, lambda = 1)
  

  ###--------------------------
  # Step4: embedding metric hn, inverse of the Riemannian matrix, symmetric
  ###--------------------------
  hn <- riemann_metric(Y = fn, laplacian = Ln, d = s, invert.h = invert.h) # array of N*s*s
  
  return(list(embedding=fn, 
              rmetric=hn, 
              weighted_matrix=W,
              adj_matrix=Kn,
              laplacian=Ln,
              nn2res = nn2res
              ))
}



# Function for graph Laplacian
# Input: W: N*N weight matrix for the neighborhood graph, bandwidth: bandwidth parameter
# Output: L: N*N graph Laplacian matrix
Laplacian <- function(W, bandwidth = 0.4, lambda = 1){
  
  D <- Matrix::Diagonal(x = rowSums(W)^(-lambda)) # inverse of a diagonal matrix
  W1 <- D %*% W %*% D
  D1 <- Matrix::Diagonal(x = 1 / rowSums(W1)) # inverse of Tn1
  L <- 4 / (bandwidth^2) * (D1 %*% W1 - Matrix::Diagonal(nrow(W))) # c=1/4 for heat kernel, depending on the choice of weights
  
  return(L)
}



# Function for Riemannian metric for each point
# The Riemannian metric and its dual associated with an embedding Y. 
# The intrinsic dimension d
# The Riemannian metric is currently denoted by G, its dual metric by H, and the Laplacian by L. 
# G at each point is the matrix inverse of H.
# `invert.h` controls whether the dual metric should be returned or not.
riemann_metric <- function(Y, laplacian, d, invert.h = FALSE){
  
  # TODO: add dimension check for all inputs
  
  H <- array(NA, dim = c(d, d, nrow(Y)))
  
  for (i in 1:d) {
    for (j in i:d) {
      yij <- Y[, i] * Y[, j]
      H[i, j, ] <- as.vector(0.5 * (laplacian %*% yij - Y[, i] * (laplacian %*% Y[, j]) - Y[, j] * (laplacian %*% Y[, i])))
      H[j, i, ] <- H[i, j, ] # symmetric matrix
    }
  }
  
  ## Pseudo inverse of H gives the final embedding metric h
  ## Array H corresponds to \tilde{H} in Step 4(a)
  ## The embedding metric H is the pseudo inverse of \tilde{H}
  for (i in 1:nrow(Y)) {
    Hsvals <- tryCatch(eigen(H[ , ,i])) # could use ginv()
    # if(class(Hsvals) != "list") browser()
    if(class(Hsvals) != "eigen") stop("Please choose a larger radius for nearest neighbor search for a connected neighborhood graph!")
    Huu <- Hsvals$vectors
    Hvv <- Hsvals$values[1:d] # top d largest eigenvalues, already sorted in decreasing order
    Hvv1 <- diag(x = 1 / Hvv, nrow = d)
    H[ , ,i] <- Huu %*% Hvv1 %*% t(Huu)
    H[, , i] <- 0.5 * (H[, , i] + t(H[, , i])) # fix not symmetric issue
  }
  
  if(invert.h){
    for (i in 1:nrow(Y)) {
      H[, , i] <- solve(H[, , i])
      H[, , i] <- 0.5 * (H[, , i] + t(H[, , i])) # fix not symmetric issue
    }
  }
  
  return(H)
}
