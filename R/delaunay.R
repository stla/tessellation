#' @title Delaunay triangulation
#' @description Delaunay triangulation (or tesellation) of a set of points.
#'
#' @param points the points given as a matrix, one point per row
#' @param atinfinity Boolean, whether to include a point at infinity
#' @param degenerate Boolean, whether to include degenerate tiles
#'
#' @return The Delaunay tessellation with many details, in a list.
#' @export
#' @useDynLib tessellation, .registration = TRUE
#' @importFrom hash hash
#' @seealso \code{\link{getDelaunaySimplicies}}
#' @examples library(tesselation)
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

#' @title Delaunay simplicies
#' @description Get Delaunay simplicies.
#'
#' @param tessellation the output of \code{\link{delaunay}}
#' @param hashes Boolean, whether to return the simplicies as hash maps
#'
#' @return The list of simplicies of the Delaunay tessellation.
#' @export
#' @importFrom hash values
#'
#' @examples library(tessellation)
#' pts <- rbind(
#'   c(-5, -5,  16),
#'   c(-5,  8,   3),
#'   c(4,  -1,   3),
#'   c(4,  -5,   7),
#'   c(4,  -1, -10),
#'   c(4,  -5, -10),
#'   c(-5,  8, -10),
#'   c(-5, -5, -10)
#' )
#'
#' tess <- delaunay(pts)
#' getDelaunaySimplicies(tess)
getDelaunaySimplicies <- function(tessellation, hashes = FALSE){
  simplicies <-
    lapply(lapply(tessellation[["tiles"]], `[[`, "simplex"), `[[`, "vertices")
  if(!hashes){
    simplicies <- lapply(simplicies, function(simplex) t(values(simplex)))
  }
  simplicies
}
