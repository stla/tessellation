#' @title Centric cuboctahedron
#' @description A cuboctahedron (12 vertices) with a point added at its center.
#'
#' @return A numeric matrix with 13 rows and 3 columns.
#' @export
centricCuboctahedron <- function(){
  rbind(
    c(-1, -1, 0),
    c(-1, 1, 0),
    c(1, -1, 0),
    c(1, 1, 0),
    c(-1, 0, -1),
    c(-1, 0, 1),
    c(1, 0, -1),
    c(1, 0, 1),
    c(0, -1, -1),
    c(0, -1, 1),
    c(0, 1, -1),
    c(0, 1, 1),
    c(0, 0, 0)
  )
}
