#' @title Sunflower
#' @description Returns a cloud of points looking like a sunflower.
#'
#' @param nmin,nmax range of the index (try on the example)
#'
#' @return A numeric matrix of 2D points (rows).
#' @export
#'
#' @examples
#' pts <- sunflower(50L, 150L)
#' opar <- par(mar = c(0, 0, 0, 0))
#' plot(pts, pch = 19, asp = 1, axes = FALSE, xlab = NA, ylab = NA)
sunflower <- function(nmin, nmax){
  n <- c(0L, nmin:nmax)
  r <- n * n
  theta <- n * pi*(3 - sqrt(5))
  x <- r * cos(theta)
  y <- r * sin(theta)
  return(cbind(x, y))
}
