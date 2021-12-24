#' @title Sunflower
#' @description Returns a cloud of points looking like a sunflower.
#'
#' @param nmax maximum value of the index (try it on the example)
#'
#' @return A numeric matrix with \code{nmax} rows and two columns.
#' @export
#'
#' @examples
#' pts <- sunflower(150L)
#' opar <- par(mar = c(0, 0, 0, 0))
#' plot(pts, pch = 19, asp = 1, axes = FALSE, xlab = NA, ylab = NA)
sunflower <- function(nmax){
  n <- 1L:nmax
  r <- sqrt(n) / 10
  theta <- n * pi*(3 - sqrt(5))
  x <- r * cos(theta)
  y <- r * sin(theta)
  return(cbind(x, y))
}
