#' t-Distributed Stochastic Neighborhood Embedding
#'
#' An S4 Class for t-SNE.
#'
#' t-SNE is a method that uses Kullback-Leibler divergence between the
#' distance matrices in high and low-dimensional space to embed the
#' data. The method is very well suited to visualize complex
#' structures in low dimensions.
#'
#' @template dimRedMethodSlots
#'
#' @template dimRedMethodGeneralUsage
#'
#' @section Parameters:
#' t-SNE can take the following parameters:
#' \describe{
#'   \item{d}{A distance function, defaults to euclidean distances}
#'   \item{perplexity}{The perplexity parameter, roughly equivalent to neighborhood size.}
#'   \item{theta}{Approximation for the nearest neighbour search, large values are more inaccurate.}
#'   \item{ndim}{The number of embedding dimensions.}
#' }
#'
#' @section Implementation:
#'
#' Wraps around \code{\link[Rtsne]{Rtsne}}, which is very well
#' documented. Setting \code{theta = 0} does a normal t-SNE, larger
#' values for \code{theta < 1} use the Barnes-Hut algorithm which
#' scales much nicer with data size. Larger values for perplexity take
#' larger neighborhoods into account.
#'
#' @references
#' Maaten, L. van der, 2014. Accelerating t-SNE using Tree-Based
#' Algorithms. Journal of Machine Learning Research 15, 3221-3245.
#'
#' van der Maaten, L., Hinton, G., 2008. Visualizing Data using
#' t-SNE. J. Mach. Learn. Res. 9, 2579-2605.
#'
#' @examples
#' \dontrun{
#' dat <- loadDataSet("3D S Curve", n = 300)
#' emb <- embed(dat, "tSNE", perplexity = 80)
#' plot(emb, type = "2vars")
#' }
#' @include dimRedResult-class.R
#' @include dimRedMethod-class.R
#' @family dimensionality reduction methods
#' @export anntSNE
#' @exportClass anntSNE
anntSNE <- setClass(
  "anntSNE",
  contains = "dimRedMethod",
  prototype = list(
    stdpars = list(nn.idx = NULL,
                   nn.dists = NULL,
                   knn = 30,
                   perplexity = 30, 
                   theta = 0.5,
                   ndim = 2,
                   annmethod = "kdtree", 
                   eps = 0,
                   radius = 1,
                   nt = 50, 
                   nlinks = 16, 
                   ef.construction = 200,
                   distance = c("euclidean", "manhattan"),
                   treetype = c("kd", "bd"), 
                   searchtype = c("standard", "priority", "radius")),
    fun = function (data, pars,
                    keep.org.data = TRUE) {
      # chckpkg("Rtsne")
      
      meta <- data@meta
      orgdata <- if (keep.org.data) data@data else NULL
      indata <- data@data
      
      if (is.null(pars$eps))      pars$eps <- 0
      # if (is.null(pars$get_geod)) pars$get_geod <- FALSE
      if (is.null(pars$distance)) pars$distance <- "euclidean"
      if (length(pars$distance) > 1) pars$distance <- pars$distance[1]
      
      if (!is.null(pars$nn.idx) & !is.null(pars$nn.dists)) {
        nn2res <- list(nn.idx = pars$nn.idx, nn.dists = pars$nn.dists)
      } else {
        nn2res <- dplyr::case_when(pars$distance=="euclidean" ~ 
                                     RANN::nn2(data = indata, query = indata, k = pars$knn + 1, 
                                               treetype = pars$treetype, searchtype = pars$searchtype, eps = pars$eps,
                                               radius = pars$radius),
                                   pars$distance=="manhattan" ~ 
                                     RANN.L1::nn2(data = indata, query = indata, k = pars$knn + 1, 
                                                  treetype = pars$treetype, searchtype = pars$searchtype, eps = pars$eps,
                                                  radius = pars$radius),
        )
        names(nn2res) <- c("nn.idx", "nn.dists")
      }
      
      
      Y_in <- matrix(runif(nrow(indata)*2), ncol = pars$ndim) # See note 2
      
      outdata <- Rtsne:::Rtsne_nn_cpp(t(nn2res$nn.index[,-1] - 1L), t(nn2res$nn.dist[,-1]), 
                           # origD=ncol(iris_matrix), 
                           no_dims = pars$ndim,
                           Y_in = Y_in,
                           init = TRUE, 
                           perplexity = pars$perplexity,
                           theta = pars$theta, 
                           max_iter = 1000,
                           verbose = TRUE,
                           stop_lying_iter = 250L, mom_switch_iter = 250L, 
                           momentum = 0.5, final_momentum = 0.8, eta = 200, exaggeration_factor = 12, num_threads = 1)$Y %>% t()
      
      
      # outdata <- Rtsne::Rtsne_neighbors(index = t(nn2res$nn.idx[,-1] - 1), 
      #                                   distance = t(nn2res$nn.dists[,-1]), 
      #                                   dims = pars$ndim, 
      #                                   perplexity = pars$perplexity,
      #                                   theta = pars$theta)$Y
      
      # outdata <- Rtsne::Rtsne(pars$d(indata),
      #                         perplexity = pars$perplexity,
      #                         theta = pars$theta,
      #                         dims = pars$ndim)$Y
      
      colnames(outdata) <- paste0("tSNE", 1:ncol(outdata))
      
      return(new(
        "dimRedResult",
        data         = new("dimRedData",
                           data = outdata,
                           meta = meta),
        org.data     = orgdata,
        has.org.data = keep.org.data,
        method       = "tsne",
        pars         = pars,
        running.time = c(0, 0, 0)
      ))
    })
)
