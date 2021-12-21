#' @title Delaunay triangulation
#' @description Delaunay triangulation (or tesellation) of a set of points.
#'
#' @param points the points given as a matrix, one point per row
#' @param atinfinity Boolean, whether to include a point at infinity
#' @param degenerate Boolean, whether to include degenerate tiles
#'
#' @return The Delaunay tessellation with many details, in a list. This list
#'   contains three fields:
#' \describe{
#'   \item{\emph{vertices}}{the vertices (or sites) of the tessellation; these are the
#'   points passed to the function}
#'   \item{\emph{tiles}}{the tiles of the tessellation (triangles in dimension 2,
#'   tetrahedra in dimension 3)}
#'   \item{\emph{tilefacets}}{the facets of the tiles of the tesselation}
#' }
#' The \strong{vertices} field is a list with the following fields:
#' \describe{
#'   \item{\emph{id}}{the id of the vertex; this is nothing but the index of
#'   the corresponding point passed to the function}
#'   \item{\emph{neighvertices}}{the ids of the vertices of the tessellation
#'   connected to this vertex by an edge}
#'   \item{\emph{neightilefacets}}{the ids of the tile facets this vertex
#'   belongs to}
#'   \item{\emph{neightiles}}{the ids of the tiles this vertex belongs to}
#' }
#' The \strong{tiles} field is a list with the following fields:
#' \describe{
#'   \item{\emph{id}}{the id of the tile}
#'   \item{\emph{simplex}}{a list describing the simplex (that is, the tile);
#'   this list contains four fields: \emph{vertices}, a
#'   \code{\link[hash]{hash}} giving the simplex vertices and their id,
#'   \emph{circumcenter}, the circumcenter of the simplex, \emph{circumradius},
#'   the circumradius of the simplex, and \emph{volume}, the volume of the
#'   simplex}
#'   \item{\emph{facets}}{the ids of the facets of this tile}
#'   \item{\emph{neighbors}}{the ids of the tiles adjacent to this tile}
#'   \item{\emph{family}}{xxxxxxx}
#'   \item{\emph{orientation}}{\code{1} or \code{-1}, an indicator of the
#'   orientation of the tile}
#' }
#' The \strong{tilefacets} field is a list with the following fields:
#' \describe{
#'   \item{\emph{id}}{the id of this tile facet}
#'   \item{\emph{subsimplex}}{a list describing the subsimplex (that is, the
#'   tile facet); this list is similar to the \emph{simplex} list of
#'   \strong{tiles}}
#'   \item{\emph{facetOf}}{one or two ids, the id(s) of the tile this facet
#'   belongs to}
#'   \item{\emph{normal}}{a vector, the normal of the tile facet}
#'   \item{\emph{offset}}{a number, the offset of the tile facet}
#' }
#'
#' @export
#' @useDynLib tessellation, .registration = TRUE
#' @importFrom hash hash
#' @seealso \code{\link{getDelaunaySimplicies}}
#' @examples library(tessellation)
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
#' del <- delaunay(points)
#' del$vertices[[1]]
#' del$tiles[[1]]
#' del$tilefacets[[1]]
delaunay <- function(points, atinfinity = FALSE, degenerate = FALSE){
  if(!is.matrix(points) || !is.numeric(points)){
    stop("The `points` argument must be a numeric matrix.", call. = TRUE)
  }
  if(ncol(points) < 2L){
    stop("The dimension must be at least 2.", call. = TRUE)
  }
  if(nrow(points) <= ncol(points)){
    stop("Insufficient number of points.", call. = TRUE)
  }
  if(any(is.na(points))){
    stop("Points with missing values are not allowed.", call. = TRUE)
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
  pointsAsList <- lapply(1L:nrow(points), function(i) points[i, ])
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
  attr(tess, "points") <- points
  tess
}

#' @title Delaunay simplicies
#' @description Get Delaunay simplicies (tiles).
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

#' @title Plot 3D Delaunay tessellation
#' @description Plot a 3D Delaunay tessellation with \strong{rgl}.
#'
#' @param tesselation the output of \code{\link{delaunay}}
#' @param color Boolean, whether to use colors
#' @param hue,luminosity if \code{color = TRUE}, these arguments are passed to
#'   \code{\link[randomcoloR]{randomColor}}
#' @param alpha opacity, number between 0 and 1
#'
#' @return No value, just renders a 3D plot.
#' @export
#' @importFrom randomcoloR randomColor
#' @importFrom utils combn
#' @importFrom rgl triangles3d
#' @importFrom hash keys values
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
#' tess <- delaunay(pts)
#' library(rgl)
#' open3d(windowRect = c(50, 50, 562, 562))
#' plotDelaunay3D(tess)
plotDelaunay3D <- function(
  tesselation, color = TRUE, hue = "random", luminosity = "light", alpha = 0.3
){
  vertices <- attr(tesselation, "points")
  if(ncol(vertices) != 3L){
    stop(
      sprintf("Invalid dimension (%d instead of 3).", ncol(vertices)),
      call. = TRUE
    )
  }
  simplicies <- getDelaunaySimplicies(tesselation, hashes = TRUE)
  edges <- unique(do.call(rbind, lapply(simplicies, function(simplex){
    t(combn(as.integer(keys(simplex)), 2L))
  })))
  nsimplicies <- length(simplicies)
  if(color){
    colors <- randomColor(nsimplicies, hue = hue, luminosity = luminosity)
    triangles <- combn(4L, 3L)
    for(i in 1L:nsimplicies){
      simplex <- t(values(simplicies[[i]]))
      for(j in 1L:4L){
        triangles3d(simplex[triangles[, j], ], color = colors[i], alpha = alpha)
      }
    }
  }
  for(i in 1L:nrow(edges)){
    edge <- edges[i, ]
    p1 <- vertices[edge[1L], ]
    p2 <- vertices[edge[2L], ]
    lines3d(rbind(p1, p2), color = "black")
  }
}

#' tile facets a vertex belongs to
#' @noRd
vertexNeighborFacets <- function(tessellation, vertexId){
  vertex <- tessellation[["vertices"]][[vertexId]]
  neighs <- vertex[["neightilefacets"]]
  tessellation[["tilefacets"]][neighs]
}
