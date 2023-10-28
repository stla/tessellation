#' @importFrom cxhull cxhull
#' @importFrom hash keys
#' @noRd
exteriorDelaunayEdges <- function(tessellation){
  tilefacets <- tessellation[["tilefacets"]]
  points <- attr(tessellation, "points")
  exteriorFacets <- Filter(Negate(sandwichedFacet), tilefacets)
  edges <- NULL
  prec <- sqrt(.Machine[["double.eps"]])
  for(tilefacet in exteriorFacets){
    normal <- tilefacet[["normal"]]
    offset <- tilefacet[["offset"]]
    # center <- tilefacet[["subsimplex"]][["circumcenter"]]
    # if(abs(crossprod(center, normal) + offset) >= sqrt(.Machine$double.eps)){
    #   normal <- -normal
    #   cat("orientation:\n")
    #
    #   tileparent <- tess$tiles[[tilefacet[["facetOf"]]]]
    #
    #   print(tileparent[["orientation"]])
    # }else{
    #   cat("nochaneg\n")
    # }
    ids <- sort(as.integer(keys(tilefacet[["subsimplex"]][["vertices"]])))
    pt1 <- points[ids[1L], ]
    pt2 <- points[ids[2L], ]
    if(
      abs(crossprod(pt1, normal) + offset) < prec
      && abs(crossprod(pt2, normal) + offset) < prec
    ){
      edges <- rbind(
        edges,
        c(ids[1L], ids[2L]),
        c(ids[2L], ids[3L]),
        c(ids[1L], ids[3L])
      )
    }
  }
  edges <- unique(edges)
  hullfacets <- cxhull(points)[["facets"]]
  edges3 <- list()
  vertices <- matrix(nrow = 0L, ncol = 3L)
  rnames <- integer(0L)
  for(i in 1L:nrow(edges)){
    edge <- edges[i, ]
    idA <- edge[1L]
    A <- points[idA, ]
    x <- vapply(hullfacets, function(f){
      c(crossprod(f[["normal"]], A)) + f[["offset"]]
    }, numeric(1L))
    Abelongs <- which(abs(x) < prec)
    idB <- edge[2L]
    B <- points[idB, ]
    if(length(Abelongs)){
      x <- vapply(Abelongs, function(j){
        f <- hullfacets[[j]]
        c(crossprod(f[["normal"]], B)) + f[["offset"]]
      }, numeric(1L))
      if(any(abs(x) < prec)){
        edges3 <-
          append(edges3, list(Edge3$new(A = A, B = B, idA = idA, idB = idB)))
        vertices <- rbind(vertices, A, B)
        rnames <- c(rnames, c(idA, idB))
      }
    }
  }
  # edges <- unique(do.call(rbind, lapply(exteriorFacets, function(f){
  #   ids <- sort(as.integer(keys(f[["subsimplex"]][["vertices"]])))
  #   rbind(
  #     c(ids[1L], ids[2L]),
  #     c(ids[2L], ids[3L]),
  #     c(ids[1L], ids[3L])
  #   )
  # })))
  # A_Bs <- apply(edges, 1L, function(edge){
  #   paste0(edge[1L], "-", edge[2L])
  # })
  # print(A_Bs)
  # unique_edges <- edges[which(table(A_Bs) == 1L), ]
  # print(unique_edges)
  #  edges < unique(edges)
  o <- order(rnames)
  vertices <- vertices[o, ]
  rownames(vertices) <- as.character(rnames[o])
  # vertices <- points[unique(c(edges)), ]
  # nedges <- nrow(edges)
  # edges3 <- vector("list", length = nedges)
  # for(i in 1L:nedges){
  #   edge <- edges[i,]
  #   edges3[[i]] <-
  #     Edge3$new(A = points[edge[1L], ], B = points[edge[2L], ])
  # }
  attr(edges3, "vertices") <- unique(vertices)
  edges3
}

volume_under_triangle <- function(x, y, z){
  sum(z) *
    (x[1L]*y[2L] - x[2L]*y[1L] + x[2L]*y[3L] -
       x[3L]*y[2L] + x[3L]*y[1L] - x[1L]*y[3L]) / 6
}

#' @title Delaunay triangulation
#' @description Delaunay triangulation (or tessellation) of a set of points.
#'
#' @param points the points given as a matrix, one point per row
#' @param atinfinity Boolean, whether to include a point at infinity;
#'   ignored if \code{elevation=TRUE}
#' @param degenerate Boolean, whether to include degenerate tiles;
#'   ignored if \code{elevation=TRUE}
#' @param exteriorEdges Boolean, for dimension 3 only, whether to return
#'   the exterior edges (see below)
#' @param elevation Boolean, only for three-dimensional points; if \code{TRUE},
#'   the function performs an elevated Delaunay triangulation (also called 2.5D
#'   Delaunay triangulation), using the third coordinate of a point as its
#'   elevation; see the example
#'
#' @return If the function performs an elevated Delaunay tessellation, then
#'   the returned value is a list with four fields: \code{mesh}, \code{edges},
#'   \code{volume}, and \code{surface}. The \code{mesh} field is an object of
#'   class \code{mesh3d}, ready for plotting with the \strong{rgl} package. The
#'   \code{edges} field is an integer matrix which provides the indices of the
#'   vertices of the edges, and an indicator of whether an edge is a border
#'   edge; this matrix is obtained with \code{\link[Rvcg]{vcgGetEdge}}.
#'   The \code{volume} field provides the sum of the
#'   volumes under the Delaunay triangles, that is to say the total volume
#'   under the triangulated surface. Finally, the \code{surface} field provides
#'   the sum of the areas of the Delaunay triangles, thus this is an approximate
#'   value of the area of the surface that is triangulated.
#'   The elevated Delaunay tessellation is built with the help of the
#'   \strong{interp} package.
#'
#' Otherwise, the function returns the Delaunay tessellation with many details,
#'   in a list. This list contains five fields:
#' \describe{
#'   \item{\emph{vertices}}{the vertices (or sites) of the tessellation; these
#'   are the points passed to the function}
#'   \item{\emph{tiles}}{the tiles of the tessellation (triangles in dimension
#'   2, tetrahedra in dimension 3)}
#'   \item{\emph{tilefacets}}{the facets of the tiles of the tessellation}
#'   \item{\emph{mesh}}{a 'rgl' mesh (\code{\link[rgl]{mesh3d}} object)}
#'   \item{\emph{edges}}{a two-columns integer matrix representing the edges,
#'   each row represents an edge; the two integers of a row are the indices of
#'   the two points which form the edge.}
#' }
#' In dimension 3, the list contains an additional field \emph{exteriorEdges}
#'   if you set \code{exteriorEdges = TRUE}. This is the list of the exterior
#'   edges, represented as \code{\link{Edge3}} objects. This field is involved
#'   in the function \code{\link{plotDelaunay3D}}.
#'
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
#'   \item{\emph{family}}{two tiles have the same family if they share the
#'   same circumcenter; in this case the family is an integer, and the family is
#'   \code{NA} for tiles which do not share their circumcenter with any other
#'   tile}
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
#' @importFrom hash hash keys
#' @importFrom rgl tmesh3d
#' @importFrom Rvcg vcgGetEdge
#'
#' @note The package provides the functions \code{\link{plotDelaunay2D}} to
#'   plot a 2D Delaunay tessellation and \code{\link{plotDelaunay3D}} to
#'   plot a 3D Delaunay tessellation. But there is no function to plot an
#'   elevated Delaunay tessellation; the examples show how to plot such a
#'   Delaunay tessellation.
#'
#' @seealso \code{\link{getDelaunaySimplices}}
#'
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
#'
#' # an elevated Delaunay tessellation ####
#' f <- function(x, y){
#'   dnorm(x) * dnorm(y)
#' }
#' x <- y <- seq(-5, 5, length.out = 50)
#' grd <- expand.grid(x = x, y = y) # grid on the xy-plane
#' points <- as.matrix(transform( # data (x_i, y_i, z_i)
#'   grd, z = f(x, y)
#' ))
#' del <- delaunay(points, elevation = TRUE)
#' del[["volume"]] # close to 1, as expected
#' # plotting
#' library(rgl)
#' mesh <- del[["mesh"]]
#' open3d(windowRect = c(100, 100, 612, 356), zoom = 0.6)
#' aspect3d(1, 1, 20)
#' shade3d(mesh, color = "limegreen", polygon_offset = 1)
#' wire3d(mesh)
#'
#' # another elevated Delaunay triangulation, to check the correctness
#' #   of the calculated surface and the calculated volume ####
#' library(Rvcg)
#' library(rgl)
#' cap <- vcgSphericalCap(angleRad = pi/2, subdivision = 3)
#' open3d(windowRect = c(100, 100, 612, 356), zoom = 0.6)
#' shade3d(cap, color = "lawngreen", polygon_offset = 1)
#' wire3d(cap)
#' # exact value of the surface of the spherical cap:
#' R <- 1
#' h <- R * (1 - sin(pi/2/2))
#' 2 * pi * R * h
#' # our approximation:
#' points <- t(cap$vb[-4, ]) # the points on the spherical cap
#' del <- delaunay(points, elevation = TRUE)
#' del[["surface"]]
#' # try to increase `subdivision` in `vcgSphericalCap` to get a
#' #   better approximation of the true value
#' # note that 'Rvcg' returns the same result as ours:
#' vcgArea(cap)
#' # let's check the volume as well:
#' pi * h^2 * (R - h/3) # true value
#' del[["volume"]]
#' # there's a warning with 'Rvcg':
#' tryCatch(vcgVolume(cap), warning = function(w) message(w))
#' suppressWarnings({vcgVolume(cap)})
delaunay <- function(
  points, atinfinity = FALSE, degenerate = FALSE, exteriorEdges = FALSE,
  elevation = FALSE
){
  stopifnot(isBoolean(atinfinity))
  stopifnot(isBoolean(degenerate))
  stopifnot(isBoolean(exteriorEdges))
  stopifnot(isBoolean(elevation))
  if(!is.matrix(points) || !is.numeric(points)){
    stop("The `points` argument must be a numeric matrix.", call. = TRUE)
  }
  dimension <- ncol(points)
  if(dimension < 2L){
    stop("The dimension must be at least 2.", call. = TRUE)
  }
  if(nrow(points) <= dimension){
    stop("Insufficient number of points.", call. = TRUE)
  }
  if(any(is.na(points))){
    stop("Points with missing values are not allowed.", call. = TRUE)
  }
  if(elevation){
    if(dimension != 3L){
      stop(
        "To get an elevated Delaunay tessellation (`elevation=TRUE`), ",
        "you have to provide three-dimensional points.",
        call. = TRUE
      )
    }
    del <- delaunay(
      points[, c(1L, 2L)], atinfinity = FALSE, degenerate = FALSE,
      exteriorEdges = FALSE, elevation = FALSE
    )

    # x <- points[, 1L]
    # y <- points[, 2L]
    # x <- (x - min(x)) / diff(range(x))
    # y <- (y - min(y)) / diff(range(y))
    # #xy <- points[, c(1L, 2L)]
    # o <- order(round(x+y, 6L), y-x)
    # xy <- cbind(x, y)[o, ]
    # if(anyDuplicated(xy)){
    #   stop("There are some duplicated points.", call. = TRUE)
    # }
    # points <- points[o, ]
    # # xy <- round(points[, c(1L, 2L)], 6)
    # Triangles <- triangles(tri.mesh(xy[, 1L], xy[, 2L]))

    delVertices <- del[["vertices"]]
    ids <- vapply(delVertices, `[[`, integer(1L), "id")

    #vertices <- do.call(cbind, lapply(delVertices, `[[`, "point"))

    vertices <- points[ids, ]

    Triangles <- do.call(rbind, lapply(del[["tiles"]], function(tile){
      indices <- tile[["vertices"]]
      if(tile[["orientation"]] == -1L){
        indices <- indices[c(2L, 1L, 3L)]
      }
      indices
    }))
    # Triangles <- Triangles[, 1L:3L]
    # vertices <- points
    mesh <- tmesh3d(
      vertices = t(vertices),
      indices = t(Triangles),
      homogeneous = FALSE
    )
    # edges <- t(vapply(del[["tilefacets"]], function(x){
    #   as.integer(keys(x[["subsimplex"]][["vertices"]]))
    # }, numeric(2L)))
    # volumes <- apply(Triangles, 1L, function(trgl){
    #   trgl <- vertices[trgl, ]
    #   volume_under_triangle(trgl[, 1L], trgl[, 2L], trgl[, 3L])
    # })
    # areas <- apply(Triangles, 1L, function(trgl){
    #   trgl <- vertices[trgl, ]
    #   triangleArea(trgl[1L, ], trgl[2L, ], trgl[3L, ])
    # })
    volumes_and_areas <- apply(Triangles, 1L, function(trgl){
      trgl <- vertices[trgl, ]
      c(
        volume_under_triangle(trgl[, 1L], trgl[, 2L], trgl[, 3L]),
        triangleArea(trgl[1L, ], trgl[2L, ], trgl[3L, ])
      )
    })
    out <- list(
      "mesh"    = mesh,
      "edges"   = `colnames<-`(
        as.matrix(vcgGetEdge(mesh))[, -3L], c("v1", "v2", "border")
      ),
      "volume"  = sum(volumes_and_areas[1L, ]),
      "surface" = sum(volumes_and_areas[2L, ])
    )
    attr(out, "elevation") <- TRUE
    class(out) <- "delaunay"
    return(out)
  }
  if(anyDuplicated(points)){
    stop("There are some duplicated points.", call. = TRUE)
  }
  storage.mode(points) <- "double"
  errfile <- tempfile(fileext = ".txt")
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
    tile <- tiles[[i]]
    simplex <- tile[["simplex"]]
    vertices <- simplex[["vertices"]]
    tess[["tiles"]][[i]] <- c(
      list("vertices" = vertices),
      tile
    )
    tess[["tiles"]][[i]][["simplex"]][["vertices"]] <-
      hash(as.character(vertices), pointsAsList[vertices])
  }
  tilefacets <- tess[["tilefacets"]]
  if(dimension == 3L){
    nTriangles <- length(tilefacets)
    Triangles <- matrix(NA_integer_, nrow = nTriangles, ncol = 3L)
    for(i in 1L:nTriangles){
      subsimplex <- tilefacets[[i]][["subsimplex"]]
      vertices <- Triangles[i, ] <- subsimplex[["vertices"]]
      tess[["tilefacets"]][[i]][["subsimplex"]][["vertices"]] <-
        hash(as.character(vertices), pointsAsList[vertices])
    }
  }else{
    for(i in seq_along(tilefacets)){
      subsimplex <- tilefacets[[i]][["subsimplex"]]
      vertices <- subsimplex[["vertices"]]
      tess[["tilefacets"]][[i]][["subsimplex"]][["vertices"]] <-
        hash(as.character(vertices), pointsAsList[vertices])
    }
  }
  attr(tess, "points") <- points
  if(dimension == 3L && exteriorEdges){
    tess[["exteriorEdges"]] <- exteriorDelaunayEdges(tess)
  }
  if(dimension == 2L){
    attr(tess[["tiles"]], "info") <-
      "Dimension 2. Tiles are triangles. A simplex volume is a triangle area."
    attr(tess[["tilefacets"]], "info") <- paste0(
      "Dimension 2. Tile facets are the edges of the triangles. ",
      "A subsimplex volume is nothing but the length of an edge."
    )
  }else if(dimension == 3L){
    tess[["mesh"]] <- mesh <- tmesh3d(
      vertices = t(points),
      indices = t(Triangles)
    )
    tess[["edges"]] <- `colnames<-`(
      as.matrix(vcgGetEdge(mesh))[, c(1L, 2L)],
      c("v1", "v2")
    )
    attr(tess[["tiles"]], "info") <-
      "Dimension 3. Tiles are tetrahedra."
    attr(tess[["tilefacets"]], "info") <- paste0(
      "Dimension 3. Tile facets are the triangles. ",
      "A subsimplex volume is nothing but the area of a triangle."
    )
  }
  class(tess) <- "delaunay"
  tess
}

#' @title Delaunay simplices
#' @description Get Delaunay simplices (tiles).
#'
#' @param tessellation the output of \code{\link{delaunay}}
#' @param hashes Boolean, whether to return the simplices as hash maps
#'
#' @return The list of simplices of the Delaunay tessellation.
#' @export
#' @importFrom hash values
#' @name getDelaunaySimplices
#' @rdname getDelaunaySimplices
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
#' getDelaunaySimplices(tess)
getDelaunaySimplices <- function(tessellation, hashes = FALSE){
  stopifnot(isBoolean(hashes))
  if(!inherits(tessellation, "delaunay")){
    stop(
      "The argument `tessellation` must be an output of the `delaunay` function.",
      call. = TRUE
    )
  }
  if(isTRUE(attr(tessellation, "elevation"))){
    stop(
      "This function is not conceived for elevated Delaunay tessellations.",
      call. = TRUE
    )
  }
  simplices <-
    lapply(lapply(tessellation[["tiles"]], `[[`, "simplex"), `[[`, "vertices")
  if(!hashes){
    simplices <- lapply(simplices, function(simplex) t(values(simplex)))
  }
  simplices
}

#' @export
#' @name getDelaunaySimplices
#' @rdname getDelaunaySimplices
getDelaunaySimplicies <- function(tessellation, hashes = FALSE) {
  getDelaunaySimplices(tessellation, hashes)
}


#' @title Plot 2D Delaunay tessellation
#' @description Plot a 2D Delaunay tessellation.
#'
#' @param tessellation the output of \code{\link{delaunay}}
#' @param border the color of the borders of the triangles; \code{NULL} for
#'   no borders
#' @param color controls the filling colors of the triangles, either
#'   \code{FALSE} for no color, \code{"random"} to use
#'   \code{\link[colorsGen]{randomColor}}, or \code{"distinct"} to use
#'   \code{\link[Polychrome]{createPalette}}
#' @param distinctArgs if \code{color = "distinct"}, a list of arguments
#'   passed to \code{\link[Polychrome]{createPalette}}
#' @param randomArgs if \code{color = "random"}, a list of arguments passed
#'   to \code{\link[colorsGen]{randomColor}}
#' @param lty,lwd graphical parameters
#' @param ... arguments passed to \code{\link{plot}}
#'
#' @return No value, just renders a 2D plot.
#' @export
#' @importFrom hash keys values
#' @importFrom graphics plot polygon par segments
#'
#' @examples # random points in a square
#' set.seed(314)
#' library(tessellation)
#' library(uniformly)
#' square <- rbind(
#'   c(-1, 1), c(1, 1), c(1, -1), c(-1, -1)
#' )
#' ptsin <- runif_in_cube(10L, d = 2L)
#' pts <- rbind(square, ptsin)
#' d <- delaunay(pts)
#' opar <- par(mar = c(0, 0, 0, 0))
#' plotDelaunay2D(
#'   d, xlab = NA, ylab = NA, asp = 1, color = "random",
#'   randomArgs = list(hue = "random", luminosity = "dark")
#' )
#' par(opar)
plotDelaunay2D <- function(
  tessellation, border = "black", color = "distinct",
  distinctArgs = list(seedcolors = c("#ff0000", "#00ff00", "#0000ff")),
  randomArgs = list(hue = "random", luminosity = "bright"),
  lty = par("lty"), lwd = par("lwd"), ...
){
  if(!inherits(tessellation, "delaunay")){
    stop(
      "The argument `tessellation` must be an output of the `delaunay` function.",
      call. = TRUE
    )
  }
  if(isTRUE(attr(tessellation, "elevation"))){
    stop(
      "This function is not conceived for elevated Delaunay tessellations.",
      call. = TRUE
    )
  }
  vertices <- attr(tessellation, "points")
  if(ncol(vertices) != 2L){
    stop(
      sprintf("Invalid dimension (%d instead of 2).", ncol(vertices)),
      call. = TRUE
    )
  }
  plot(vertices, type = "n", ...)
  if(!isFALSE(color)){
    color <- match.arg(color, c("random", "distinct"))
    simplices <- getDelaunaySimplices(tessellation, hashes = TRUE)
    nsimplices <- length(simplices)
    if(color == "random") {
      colors <- rcolors(nsimplices, randomArgs)
    } else {
      colors <- distinctColors(nsimplices, distinctArgs)
    }
    for(i in 1L:nsimplices){
      triangle <- t(values(simplices[[i]]))
      polygon(triangle, border = NA, col = colors[i])
    }
  }
  if(!is.null(border)){
    edges <- lapply(tessellation[["tilefacets"]], function(tilefacet){
      as.integer(keys(tilefacet[["subsimplex"]][["vertices"]]))
    })
    # cat("nedges:\n")
    # print(length(edges))
    # edges <- uniqueWith(edges, sameSegments)
    # print(length(edges))
    for(edge in edges){
      p0 <- vertices[edge[1L], ]
      p1 <- vertices[edge[2L], ]
      segments(
        p0[1L], p0[2L], p1[1L], p1[2L], col = border, lty = lty, lwd = lwd
      )
    }
  }
}


#' @title Plot 3D Delaunay tessellation
#' @description Plot a 3D Delaunay tessellation with \strong{rgl}.
#'
#' @param tessellation the output of \code{\link{delaunay}}
#' @param color controls the filling colors of the tetrahedra, either
#'   \code{FALSE} for no color, \code{"random"} to use
#'   \code{\link[colorsGen]{randomColor}}, or \code{"distinct"} to use
#'   \code{\link[Polychrome]{createPalette}}
#' @param distinctArgs if \code{color = "distinct"}, a list of arguments
#'   passed to \code{\link[Polychrome]{createPalette}}
#' @param randomArgs if \code{color = "random"}, a list of arguments passed
#'   to \code{\link[colorsGen]{randomColor}}
#' @param alpha opacity, number between 0 and 1
#' @param exteriorEdgesAsTubes Boolean, whether to plot the exterior edges
#'   as tubes; in order to use this feature, you need to set
#'   \code{exteriorEdges = TRUE} in the \code{\link{delaunay}} function
#' @param tubeRadius if \code{exteriorEdgesAsTubes = TRUE}, the radius of
#'   the tubes
#' @param tubeColor if \code{exteriorEdgesAsTubes = TRUE}, the color of
#'   the tubes
#'
#' @return No value, just renders a 3D plot.
#' @export
#' @importFrom utils combn
#' @importFrom rgl triangles3d spheres3d
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
#' plotDelaunay3D(tess, color = "random")
#' \donttest{open3d(windowRect = c(50, 50, 562, 562))
#' plotDelaunay3D(
#'   tess, exteriorEdgesAsTubes = TRUE, tubeRadius = 0.3, tubeColor = "yellow"
#' )}
plotDelaunay3D <- function(
  tessellation, color = "distinct",
  distinctArgs = list(seedcolors = c("#ff0000", "#00ff00", "#0000ff")),
  randomArgs = list(hue = "random", luminosity = "bright"),
  alpha = 0.3, exteriorEdgesAsTubes = FALSE, tubeRadius, tubeColor
){
  stopifnot(isBoolean(exteriorEdgesAsTubes))
  if(!inherits(tessellation, "delaunay")){
    stop(
      "The argument `tessellation` must be an output of the `delaunay` function.",
      call. = TRUE
    )
  }
  if(isTRUE(attr(tessellation, "elevation"))){
    stop(
      "This function is not conceived for elevated Delaunay tessellations.",
      call. = TRUE
    )
  }
  vertices <- attr(tessellation, "points")
  if(ncol(vertices) != 3L){
    stop(
      sprintf("Invalid dimension (%d instead of 3).", ncol(vertices)),
      call. = TRUE
    )
  }
  # simplices <- getDelaunaySimplices(tessellation, hashes = TRUE)
  # # edges <- unique(do.call(rbind, lapply(simplices, function(simplex){
  # #   t(combn(as.integer(keys(simplex)), 2L))
  # # })))
  # nsimplices <- length(simplices)
  if(!isFALSE(color)) {
    color <- match.arg(color, c("random", "distinct"))
    simplices <- getDelaunaySimplices(tessellation, hashes = TRUE)
    nsimplices <- length(simplices)
    if(color == "random") {
      colors <- rcolors(nsimplices, randomArgs)
    } else {
      colors <- distinctColors(nsimplices, distinctArgs)
    }
    triangles <- combn(4L, 3L)
    for(i in 1L:nsimplices) {
      simplex <- t(values(simplices[[i]]))
      for(j in 1L:4L) {
        triangles3d(simplex[triangles[, j], ], color = colors[i], alpha = alpha)
      }
    }
  }
  edges <- tessellation[["edges"]]
  for(i in 1L:nrow(edges)){
    edge <- edges[i, ]
    p1 <- vertices[edge[1L], ]
    p2 <- vertices[edge[2L], ]
    lines3d(rbind(p1, p2), color = "black")
  }
  if(exteriorEdgesAsTubes){
    if(!"exteriorEdges" %in% names(tessellation)){
      warning(
        "You didn't set the option `exteriorEdges=TRUE` in the `delaunay` ",
        "function, therefore the option `exteriorEdgesAsTubes` is ignored."
      )
    }
    edges3 <- tessellation[["exteriorEdges"]]
    for(edge3 in edges3){
      edge3$plot(
        edgeAsTube = TRUE, tubeRadius = tubeRadius, tubeColor = tubeColor
      )
    }
    spheres3d(
      attr(edges3, "vertices"), radius = 1.5*tubeRadius, color = tubeColor
    )
  }
}

#' tile facets a vertex belongs to
#' @noRd
vertexNeighborFacets <- function(tessellation, vertexId){
  vertex <- tessellation[["vertices"]][[vertexId]]
  neighs <- vertex[["neightilefacets"]]
  tessellation[["tilefacets"]][neighs]
}

#' @title Tessellation volume
#' @description The volume of the Delaunay tessellation, that is, the volume of
#'   the convex hull of the sites.
#'
#' @param tessellation output of \code{\link{delaunay}}
#'
#' @return A number, the volume of the Delaunay tessellation (area in 2D).
#' @seealso \code{\link{surface}}
#' @export
volume <- function(tessellation){
  if(!inherits(tessellation, "delaunay")){
    stop(
      "The argument `tessellation` must be an output of the `delaunay` function.",
      call. = TRUE
    )
  }
  if(isTRUE(attr(tessellation, "elevation"))){
    stop(
      "This function is not conceived for elevated Delaunay tessellations.",
      call. = TRUE
    )
  }
  tileVolumes <- vapply(tessellation[["tiles"]], function(tile){
    tile[["simplex"]][["volume"]]
  }, numeric(1L))
  sum(tileVolumes)
}

sandwichedFacet <- function(tilefacet){
  length(tilefacet[["facetOf"]]) == 2L
}

#' @title Tessellation surface
#' @description Exterior surface of the Delaunay tessellation.
#'
#' @param tessellation output of \code{\link{delaunay}}
#'
#' @return A number, the exterior surface of the Delaunay tessellation
#'   (perimeter in 2D).
#' @seealso \code{\link{volume}}
#' @note It is not guaranteed that this function provides the correct result
#'   for all cases. The exterior surface of the Delaunay tessellation is the
#'   exterior surface of the convex hull of the sites (the points), and you can
#'   get it with the \strong{cxhull} package (by summing the volumes of the
#'   facets). Moreover, I encountered some cases for which I got a correct
#'   result only with the option \code{degenerate=TRUE} in the
#'   \code{delaunay} function. I will probably remove this function in the
#'   next version.
#' @export
surface <- function(tessellation){
  if(!inherits(tessellation, "delaunay")){
    stop(
      "The argument `tessellation` must be an output of the `delaunay` function.",
      call. = TRUE
    )
  }
  if(isTRUE(attr(tessellation, "elevation"))){
    stop(
      "This function is not conceived for elevated Delaunay tessellations.",
      call. = TRUE
    )
  }
  exteriorFacets <-
    Filter(Negate(sandwichedFacet), tessellation[["tilefacets"]])
  ridgeSurfaces <- vapply(exteriorFacets, function(tilefacet){
    tilefacet[["subsimplex"]][["volume"]]
  }, numeric(1L))
  sum(ridgeSurfaces)
}
