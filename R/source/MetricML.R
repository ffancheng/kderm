#' Metric manifold learning algorithm modified from the python megaman package
#'
#' @param d A set of \code(x) data points in \mathcal{R}^r
#' @param s The number of the dimensions of the embedding
#' @param epsilon The bandwidth parameter for the manifold learning algorithm
#' @param method The manifold learning algorithm to be applied to the data \code{d}
#' @return The list of the embedding coordinates \code{fn_p} and the embedding metric \code{hn_p} for each point \code{p} \in \code{d}
#' @examples
#' add(1, 1)
#' add(10, 1)
metricML <- function(d, s, epsilon, method){
  
  # Step1: similarity matrix
  
  
  # Step2: Laplacian matrix
  
  # Step3: embedding coordinates fn_p
  
  # Step4: embedding metric hn_p
  
  return(list(fn_p, hn_p))
}
