# add argument h
# https://github.com/robjhyndman/hdrcde/blob/f767b441e4338043258bcbcf8a60b281d1fd22c9/R/hdrscatterplot.R#L43
hdrscatterplot_new <- function(x, y, levels = c(1, 50, 99), kde.package = c("ash", "ks"), noutliers = NULL, label = NULL, h = NULL, den = NULL, ...) {
  levels <- sort(levels)
  if (missing(y)) {
    data <- x
  } else {
    data <- data.frame(x = x, y = y)
    names(data) <- make.names(c(deparse(substitute(x)), deparse(substitute(y))))
  }
  
  vnames <- names(data)
  
  if(is.null(den)){
    den <- hdr.2d(data[, 1], data[, 2], prob = levels, kde.package = kde.package, h = h, ...) # add all arguments, h for bandwidth
  }
  
  region <- numeric(NROW(data)) + 100
  for (i in seq_along(levels))
    region[den$fxy > den$falpha[i]] <- 100 - levels[i]
  if (is.null(noutliers)) {
    noutliers <- sum(region > max(levels))
  }
  
  noutliers <- min(noutliers, NROW(data))
  
  
  xlim <- diff(range(data[, 1]))
  ylim <- diff(range(data[, 2]))
  
  # Construct region factor
  levels <- sort(unique(region[region < 100]), decreasing = TRUE)
  levels <- c(levels, 100)
  data$Region <- factor(region,
                        levels = levels,
                        labels = c(paste(head(levels, -1)), ">99")
  )
  # Sort data so the larger regions go first (other than outliers)
  k <- region
  k[region == 100] <- 0
  ord <- order(k, decreasing = TRUE)
  data <- data[ord, ] # data is reordered by ord
  
  if (noutliers > 0) {
    outlier_rank <- order(den$fxy[ord])  # order den$fxy as well to match the new data
    outliers <- outlier_rank[seq_len(noutliers)] # take top noutliers labels for annotation
  }
  
  p <- ggplot2::ggplot(data, ggplot2::aes_string(vnames[1], vnames[2])) +
    ggplot2::geom_point(ggplot2::aes_string(col = "Region"))
  p <- p + ggplot2::scale_colour_manual(
    name = "HDRs",
    breaks = c(paste(head(sort(levels), -1)), ">99"),
    values = c(RColorBrewer::brewer.pal(length(levels), "YlOrRd")[-1], "#000000")
  )
  
  # Show outliers
  if (is.null(label)) {
    label <- rownames(data)[outliers]
  } else {
    if (length(label) != nrow(data)) stop("The length of label is not the same as x and y!")
    label <- label[ord[outliers]]
  }
  if (noutliers > 0) {
    p <- p + ggplot2::annotate("text",
                               x = data[outliers, 1] + xlim / 50, y = data[outliers, 2] + ylim / 50,
                               label = label, col = "blue", cex = 2.5
    )
  }
  # return(p)
  return(list(p=p, 
              # outlier= rownames(data)[outlier_rank], densities = den$fxy[ord]
              outlier = order(den$fxy, decreasing = FALSE),  # densities ascending, top anomalous to typical
              densities = den$fxy, # density order matching the data input
              den = den # list of x, y, z(matrix)
              )
         ) ## TODO: check outlier/density order
}


# Main function. Called by hdr.boxplot.2d
#' @rdname hdr.boxplot.2d
#' @param den Bivariate density estimate (a list with elements x, y and z where
#' x and y are grid values and z is a matrix of density values). If
#' \code{NULL}, the density is estimated.
#' @export
hdr.2d <- function(x, y, prob = c(50, 95, 99), den=NULL, kde.package=c("ash","ks"), h=NULL,
                   xextend=0.15, yextend=0.15)
{
  # Convert prob to coverage percentage if necessary
  if(max(prob) > 50) # Assume prob is coverage percentage
    alpha <- (100-prob)/100
  else # prob is tail probability (for backward compatibility)
    alpha <- prob
  alpha <- sort(alpha)
  
  # Estimate bivariate density
  if(is.null(den))
    den <- hdrcde:::den.estimate.2d(x,y,kde.package,h,xextend,yextend)
  
  # Calculates falpha needed to compute HDR of bivariate density den.
  # Also finds approximate mode.
  fxy <- hdrcde:::interp.2d(den$x,den$y,den$z,x,y)
  falpha <- quantile(fxy, alpha)
  index <- which.max(fxy)
  mode <- c(x[index],y[index])
  return(structure(list(mode=mode,falpha=falpha,fxy=fxy, den=den, alpha=alpha, x=x, y=y), class="hdr2d"))
}


# Bivariate density estimate
# Modified from hdrcde package https://github.com/robjhyndman/hdrcde/blob/master/R/hdr.boxplot.2d.R

den.estimate.2d <- function(x, y, kde.package=c("ash","ks"), h=NULL, xextend=0.15,yextend=0.15) {
  kde.package <- match.arg(kde.package)
  # Find ranges for estimates
  xr <- diff(range(x,na.rm=TRUE))
  yr <- diff(range(y,na.rm=TRUE))
  xr <- c(min(x)-xr*xextend,max(x)+xr*xextend)
  yr <- c(min(y)-yr*yextend,max(y)+yr*yextend)
  if(kde.package=="ash")
  {
    if(is.null(h))
      h <- c(5,5)
    den <- ash::ash2(ash::bin2(cbind(x,y),rbind(xr,yr)),h)
  } else
  {
    X <- cbind(x,y)
    if(is.null(h))
      h <- ks::Hpi.diag(X,binned=TRUE) # 2*2 diagonal matrix
    else if(is.vector(h)) # input h as a vector, which was default
      h <- diag(h) 
    else
      # h <- diag(h) 
      # den <- ks::kde(x=X,H=h,xmin=c(xr[1],yr[1]),xmax=c(xr[2],yr[2])) # h is not variable
      den <- ks:::kde.sp.2d(x=X,H=h,xmin=c(xr[1],yr[1]),xmax=c(xr[2],yr[2])) 
    den <- list(x=den$eval.points[[1]],y=den$eval.points[[2]],z=den$estimate)
  }
  return(den)
}
