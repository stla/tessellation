#' Title
#'
#' @param points
#' @param atinfinity
#' @param degenerate
#'
#' @return
#' @export
#' @useDynLib tessellation, .registration = TRUE
#' @importFrom hash hash
#' @examples
#' points <- rbind(
#'  c(0.5,0.5,0.5),
#'  c(0,0,0),
#'  c(0,0,1),
#'  c(0,1,0),
#'  c(0,1,1),
#'  c(1,0,0),
#'  c(1,0,1),
#'  c(1,1,0),
#'  c(1,1,1)
#' )
#' delaunay(points)
delaunay <- function(points, atinfinity = FALSE, degenerate = FALSE){
  if(!is.matrix(points) || !is.numeric(points)){
    stop("`points` must be a numeric matrix")
  }
  if(ncol(points) < 2L){
    stop("dimension must be at least 2")
  }
  if(nrow(points) <= ncol(points)){
    stop("insufficient number of points")
  }
  if(any(is.na(points))){
    stop("missing values are not allowed")
  }
  errfile <- tempfile(fileext=".txt")
  tess <- tryCatch({
    .Call(
      "delaunay_",
      points,
      as.integer(atinfinity),
      as.integer(degenerate),
      0,
      errfile
    )
  }, error = function(e){
    try(cat(readLines(errfile), sep="\n"), silent = TRUE)
    stop(e)
  })
  pointsAsList <- lapply(1:nrow(points), function(i) points[i, ])
  tiles <- tess[["tiles"]]
  for(i in seq_along(tiles)){
    simplex <- tiles[[i]][["simplex"]]
    vertices <- simplex[["vertices"]]
    tess[["tiles"]][[i]][["simplex"]][["vertices"]] <-
      hash(as.character(vertices), pointsAsList[vertices])
  }
  tilefacets <- tess[["tilefacets"]]
  for(i in seq_along(tilefacets)){
    subsimplex <- tilefacets[[i]][["subsimplex"]]
    vertices <- subsimplex[["vertices"]]
    tess[["tilefacets"]][[i]][["subsimplex"]][["vertices"]] <-
      hash(as.character(vertices), pointsAsList[vertices])
  }
  tess
}
