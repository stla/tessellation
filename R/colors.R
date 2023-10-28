#' @importFrom Polychrome createPalette
#' @noRd
distinctColors <- function(n, argsList) {
  f <- function(...) {
    createPalette(n, ...)
  }
  do.call(f, argsList)
}

#' @importFrom colorsGen randomColor
rcolors <- function(n, argsList) {
  f <- function(...) {
    randomColor(n, ...)
  }
  do.call(f, argsList)
}
