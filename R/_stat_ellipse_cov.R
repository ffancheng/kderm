#' Compute data ellipses using variable covariance matrix for each point
#'
#' Each data point has an ellipse calculated based on a variable covariance matrix. The covariance matrix for each point is saved in a nested data frame from `tidyr::nest()`. 
#' The method for calculating the ellipses has been modified from
#' `car::dataEllipse` (Fox and Weisberg, 2011)
#' The code has been modified from `ggplot2::stat_ellipse()` where covariance matrix can't be specified. https://github.com/tidyverse/ggplot2/blob/master/R/stat-ellipse.R
#'
#' @references John Fox and Sanford Weisberg (2011). An \R Companion to
#'   Applied Regression, Second Edition. Thousand Oaks CA: Sage. URL:
#'   \url{https://socialsciences.mcmaster.ca/jfox/Books/Companion/}
#' @param level The level at which to draw an ellipse,
#'   or, if `type="euclid"`, the radius of the circle to be drawn.
#' @param type The type of ellipse.
#'   The default `"t"` assumes a multivariate t-distribution, and
#'   `"norm"` assumes a multivariate normal distribution.
#'   `"euclid"` draws a circle with the radius equal to `level`,
#'   representing the euclidean distance from the center.
#'   This ellipse probably won't appear circular unless `coord_fixed()` is applied.
#' @param segments The number of segments to be used in drawing the ellipse.
#' @inheritParams layer
#' @inheritParams geom_point
#' @export
#' @examples
#' ggplot(faithful, aes(waiting, eruptions)) +
#'   geom_point() +
#'   stat_ellipse()
#'
#' ggplot(faithful, aes(waiting, eruptions, color = eruptions > 3)) +
#'   geom_point() +
#'   stat_ellipse()
#'
#' ggplot(faithful, aes(waiting, eruptions, color = eruptions > 3)) +
#'   geom_point() +
#'   stat_ellipse(type = "norm", linetype = 2) +
#'   stat_ellipse(type = "t")
#'
#' ggplot(faithful, aes(waiting, eruptions, color = eruptions > 3)) +
#'   geom_point() +
#'   stat_ellipse(type = "norm", linetype = 2) +
#'   stat_ellipse(type = "euclid", level = 3) +
#'   coord_fixed()
#'
#' ggplot(faithful, aes(waiting, eruptions, fill = eruptions > 3)) +
#'   stat_ellipse(geom = "polygon")
stat_ellipse_cov <- function(mapping = NULL, data = NULL,
                         geom = "path", position = "identity",
                         ...,
                         type = "t",
                         level = 0.95,
                         segments = 51,
                         n.plot = 50,
                         na.rm = FALSE,
                         show.legend = NA,
                         inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatEllipseCov,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      type = type,
      level = level,
      segments = segments,
      n.plot = n.plot,
      na.rm = na.rm,
      ...
    )
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
StatEllipseCov <- ggproto("StatEllipseCov", Stat,
                       required_aes = c("x", "y", "covmat"),
                       
                       compute_group = function(data, scales, type = "t", level = 0.95,
                                                segments = 51, n.plot = 50, na.rm = FALSE) {
                         
                         calculate_ellipse_cov(data = data, vars = c("x", "y", "covmat"), type = type,
                                           level = level, segments = segments, n.plot)
                       }
)

calculate_ellipse_cov <- function(data, vars, type, level, segments, n.plot){
  dfn <- 2
  dfd <- nrow(data) - 1
  
  x <- data[, vars]
  # dfd <- nrow(x) - 1
  center <- x[, 1:2] #%>% cbind(NA, NA)
  shape <- matrix(x$covmat %>% unlist, ncol = 4, byrow = T) %>% as_tibble()
  # purrr::map2_dfc(center, shape, cal_single_ellipse_cov)

  ell <- NULL
  for (i in sample(1:(dfd+1), n.plot)) {
  temp <- cal_single_ellipse_cov(center = center[i, ], shape = shape[i, ], type, level, segments, dfd)
  ell <- rbind(ell, temp)
} 

  return(ell)
  
}

cal_single_ellipse_cov <- function(center, shape, type,
                                   level, segments, dfd) {
  dfn <- 2
  # dfd <- nrow(data) - 1
  
  # center <- as.numeric(x[, vars[1:2]])
  # shape <- x[, vars(3)][[1]]
  
  shape <- matrix(as.numeric(shape), ncol = 2, byrow = T)
  center <- as.numeric(center[1, 1:2])
  
  chol_decomp <- chol(shape)
  if (type == "euclid") {
    radius <- level/max(chol_decomp) / 50
  } else {
    radius <- sqrt(dfn * stats::qf(level, dfn, dfd)) / 50
  }
  angles <- (0:segments) * 2 * pi/segments
  unit.circle <- cbind(cos(angles), sin(angles))
  ellipse <- sweep(radius * (unit.circle %*% chol_decomp), 2, center, "+")
  
  colnames(ellipse) <- c("x", "y")
  ggplot2:::mat_2_df(ellipse)
}

