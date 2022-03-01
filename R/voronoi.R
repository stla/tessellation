sameFamily <- function(tile1, tile2){
  family1 <- tile1[["family"]]
  family2 <- tile2[["family"]]
  if(is.na(family1) || is.na(family2)){
    return(FALSE)
  }
  family1 == family2
}

edgeFromTileFacet <- function(tessellation, tilefacet){
  tileindices <- tilefacet[["facetOf"]]
  tiles <- tessellation[["tiles"]]
  tile1 <- tiles[[tileindices[1L]]]
  c1 <- tile1[["simplex"]][["circumcenter"]]
  if(any(is.nan(c1))){
    return(NULL)
  }
  if(length(tileindices) == 1L){
    return(newIEdge(c1, tilefacet[["normal"]]))
  }
  tile2 <- tiles[[tileindices[2L]]]
  c2 <- tile2[["simplex"]][["circumcenter"]]
  if(
    sameFamily(tile1, tile2) || isTRUE(all.equal(c1, c2)) || any(is.nan(c2))
  ){
    return(NULL)
  }
  newEdge(c1, c2)
}

voronoiCell <- function(facetsQuotienter, edgeTransformer){
  function(tessellation, vertexId){
    tilefacets <- facetsQuotienter(
      unname(vertexNeighborFacets(tessellation, vertexId))
    )
    edges <- Filter(Negate(is.null), lapply(tilefacets, function(tilefacet){
      edgeFromTileFacet(tessellation, tilefacet)
    }))
    nedges <- length(edges)
    bounded <- nedges > 0L
    while(nedges > 0L){
      bounded <- bounded && inherits(edges[[nedges]], c("Edge2", "Edge3"))
      nedges <- nedges - 1L
    }
    cell <- edgeTransformer(edges)
    attr(cell, "bounded") <- bounded
    cell
  }
}

#' @importFrom sets pair
#' @noRd
zip <- function(matrix, list){
  lapply(seq_len(nrow(matrix)), function(i){
    pair(site = matrix[i, ], cell = list[[i]])
  })
}

voronoi0 <- function(cellGetter, tessellation){
  vertices <- attr(tessellation, "points")
  bounded <- logical(nrow(vertices))
  L <- lapply(1L:nrow(vertices), function(i){
    cell <- cellGetter(tessellation, i)
    bounded[i] <<- attr(cell, "bounded")
    cell
  })
  out <- zip(vertices, L)
  attr(out, "nbounded") <- sum(bounded)
  out
}

#' @title Voronoï tessellation
#' @description Voronoï tessellation from Delaunay tessellation; this is a list
#'   of pairs made of a site (a vertex) and a list of edges.
#'
#' @param tessellation output of \code{\link{delaunay}}
#'
#' @return A list of pairs representing the Voronoï tessellation. Each
#'   \code{\link[sets]{pair}} is named: the first component is called
#'   \code{"site"}, and the second component is called \code{"cell"}.
#'
#' @seealso \code{\link{isBoundedCell}}, \code{\link{cellVertices}},
#'   \code{\link{plotBoundedCell2D}}, \code{\link{plotBoundedCell3D}}
#' @export
#' @importFrom english english
#'
#' @examples library(tessellation)
#' d <- delaunay(centricCuboctahedron())
#' v <- voronoi(d)
#' # the Voronoï diagram has 13 cells (one for each site):
#' length(v)
#' # there is only one bounded cell:
#' length(Filter(isBoundedCell, v)) # or attr(v, "nbounded")
voronoi <- function(tessellation){
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
  v <- voronoi0(voronoiCell(identity, identity), tessellation)
  nbounded <- attr(v, "nbounded")
  message(
    sprintf(
      "Vorono\u00ef diagram with %s bounded cell%s.",
      ifelse(nbounded <= 100L, as.character(english(nbounded)), nbounded),
      ifelse(nbounded > 1L, "s", "")
    )
  )
  class(v) <- "voronoi"
  v
}

#' @title Is this cell bounded?
#' @description Check whether a Voronoï cell is bounded, i.e. contains only
#'   finite edges.
#'
#' @param cell a Voronoï cell
#'
#' @return A Boolean value, whether the cell is bounded.
#' @export
#' @importFrom sets is.tuple
isBoundedCell <- function(cell){
  if(is.tuple(cell)){
    cell <- cell[["cell"]]
  }
  length(cell) != 0L && all(vapply(cell, function(edge){
    inherits(edge, "Edge2") || inherits(edge, "Edge3")
  }, logical(1L)))
}

#' @title Vertices of a bounded cell
#' @description Get all vertices of a bounded cell, without duplicates.
#'
#' @param cell a bounded Voronoï cell
#' @param check.bounded Boolean, whether to check that the cell is bounded;
#'   set to \code{FALSE} for a small speed gain if you know that the cell is
#'   bounded
#'
#' @return A matrix, each row represents a vertex.
#' @export
#' @importFrom sets is.tuple
#'
#' @examples library(tessellation)
#' d <- delaunay(centricCuboctahedron())
#' v <- voronoi(d)
#' cell13 <- v[[13]]
#' isBoundedCell(cell13) # TRUE
#' library(rgl)
#' open3d(windowRect = c(50, 50, 562, 562))
#' invisible(lapply(cell13[["cell"]], function(edge){
#'   edge$plot(edgeAsTube = TRUE, tubeRadius = 0.025, tubeColor = "yellow")
#' }))
#' cellvertices <- cellVertices(cell13)
#' spheres3d(cellvertices, radius = 0.1, color = "green")
cellVertices <- function(cell, check.bounded = TRUE){
  stopifnot(isBoolean(check.bounded))
  if(check.bounded && !isBoundedCell(cell)){
    stop(
      "This function applies to bounded cells only.",
      call. = TRUE
    )
  }
  if(is.tuple(cell)){
    cell <- cell[["cell"]]
  }
  stackedVertices <- lapply(cell, function(edge) edge$stack())
  vertices <- unique(do.call(rbind, stackedVertices))
  if(ncol(vertices) == 2L){ # i.e. dimension 2
    # we must order the vertices:
    vectors <- sweep(vertices, 2L, colMeans(vertices), check.margin = FALSE)
    vector1 <- vectors[1L, ]
    a <- atan2(vector1[2L], vector1[1L])
    vectors <- vectors[-1L, ]
    angles <- c(0, apply(vectors, 1L, function(v) atan2(v[2L], v[1L]) - a))
    vertices <- vertices[order(angles %% (2*pi)), ]
  }
  vertices
}

#' @title Plot a bounded Voronoï 3D cell
#' @description Plot a bounded Voronoï 3D cell with \strong{rgl}.
#'
#' @param cell a bounded Voronoï 3D cell
#' @param edgesAsTubes Boolean, whether to plot edges as tubes or as lines
#' @param tubeRadius radius of the tubes if \code{edgesAsTubes = TRUE}
#' @param tubeColor color of the tubes if \code{edgesAsTubes = TRUE}
#' @param facetsColor color of the facets; \code{NA} for no color
#' @param alpha opacity of the facets, a number between 0 and 1
#' @param check.bounded Boolean, whether to check that the cell is bounded;
#'   set to \code{FALSE} for a small speed gain if you know that the cell is
#'   bounded
#'
#' @return No value, this function just plots the cell.
#' @export
#' @importFrom cxhull cxhull
#' @importFrom rgl lines3d cylinder3d shade3d triangles3d spheres3d
#' @importFrom sets is.tuple
#'
#' @examples library(tessellation)
#' d <- delaunay(centricCuboctahedron())
#' v <- voronoi(d)
#' cell13 <- v[[13]]
#' isBoundedCell(cell13) # TRUE
#' \donttest{library(rgl)
#' open3d(windowRect = c(50, 50, 562, 562))
#' plotBoundedCell3D(
#'   cell13, edgesAsTubes = TRUE, tubeRadius = 0.03, tubeColor = "yellow",
#'   facetsColor = "navy", alpha = 0.7
#' )}
plotBoundedCell3D <- function(
  cell, edgesAsTubes = FALSE, tubeRadius, tubeColor,
  facetsColor = NA, alpha = 1, check.bounded = TRUE
){
  stopifnot(isBoolean(edgesAsTubes))
  stopifnot(isBoolean(check.bounded))
  if(check.bounded && !isBoundedCell(cell)){
    stop(
      "This function applies to bounded cells only.",
      call. = TRUE
    )
  }
  if(edgesAsTubes && (missing(tubeRadius) || missing(tubeColor))){
    stop(
      "If you use the option `edgesAsTubes`, you have to set a value for ",
      "`tubRadius` and `tubeColor`.",
      call. = TRUE
    )
  }
  if(is.tuple(cell)){
    cell <- cell[["cell"]]
  }
  if(edgesAsTubes || !is.na(facetsColor)){
    stackedVertices <- lapply(cell, function(edge) edge$stack())
    vertices <- unique(do.call(rbind, stackedVertices))
  }
  if(edgesAsTubes){
    for(edge in cell){
      tube <- cylinder3d(
        rbind(
          edge$A, edge$B
        ),
        radius = tubeRadius,
        sides = 90
      )
      shade3d(tube, color = tubeColor)
    }
    spheres3d(vertices, radius = 1.5*tubeRadius, color = tubeColor)
  }else{
    for(edge in cell){
      edge$plot()
    }
  }
  if(!is.na(facetsColor)){
    h <- cxhull(vertices, triangulate = TRUE)
    for(facet in h$facets){
      triangle <- t(vapply(
        facet$vertices, function(id) h$vertices[[as.character(id)]]$point,
        numeric(3L)
      ))
      triangles3d(triangle, color = facetsColor, alpha = alpha)
    }
  }
}


#' @title Plot a bounded Voronoï 2D cell
#' @description Plot a bounded Voronoï 2D cell.
#'
#' @param cell a bounded Voronoï 2D cell
#' @param border color of the borders of the cell; \code{NA} for no color
#' @param color color of the cell; \code{NA} for no color
#' @param check.bounded Boolean, whether to check that the cell is bounded;
#'   set to \code{FALSE} for a small speed gain if you know that the cell is
#'   bounded
#' @param ... graphical parameters for the borders
#'
#' @return No value, this function just plots the cell (more precisely, it adds
#'   the plot of the cell to the current plot).
#' @export
#' @importFrom graphics polygon segments
#' @importFrom sets is.tuple
#'
#' @examples library(tessellation)
#' centricSquare <- rbind(
#'   c(-1, 1), c(1, 1), c(1, -1), c(-1, -1), c(0, 0)
#' )
#' d <- delaunay(centricSquare)
#' v <- voronoi(d)
#' cell5 <- v[[5]]
#' isBoundedCell(cell5) # TRUE
#' plot(centricSquare, type = "n", asp = 1, xlab = "x", ylab = "y")
#' plotBoundedCell2D(cell5, color = "pink")
plotBoundedCell2D <- function(
  cell, border = "black", color = NA, check.bounded = TRUE, ...
){
  # if(!isBoundedCell(cell)){
  #   stop(
  #     "This function applies to bounded cells only.",
  #     call. = TRUE
  #   )
  # }
  # if(is.tuple(cell)){
  #   cell <- cell[["cell"]]
  # }
  stopifnot(isBoolean(check.bounded))
  if(!is.na(color)){
    vertices <- cellVertices(cell, check.bounded = check.bounded)
    # stackedVertices <- lapply(cell, function(edge) edge$stack())
    # vertices <- unique(do.call(rbind, stackedVertices))
    # # we must order the vertices:
    # vectors <- sweep(vertices, 2L, colMeans(vertices), check.margin = FALSE)
    # vector1 <- vectors[1L, ]
    # a <- atan2(vector1[2L], vector1[1L])
    # vectors <- vectors[-1L, ]
    # angles <- c(0, apply(vectors, 1L, function(v) atan2(v[2L], v[1L]) - a))
    # vertices <- vertices[order(angles %% (2*pi)), ]
    polygon(vertices, border = NA, col = color)
  }
  if(!is.na(border)){
    if(is.na(color) && check.bounded){
      if(!isBoundedCell(cell)){
        stop(
          "This function applies to bounded cells only.",
          call. = TRUE
        )
      }
    }
    if(is.tuple(cell)){
      cell <- cell[["cell"]]
    }
    for(edge in cell){
      edge$plot(color = border, ...)
    }
  }
}


#' @title Plot Voronoï diagram
#' @description Plot all the bounded cells of a 2D or 3D Voronoï tessellation.
#'
#' @param v an output of \code{\link{voronoi}}
#' @param colors this can be \code{"random"} to use random colors for the cells
#'   (with \code{\link[randomcoloR]{randomColor}}), \code{"distinct"} to use
#'   distinct colors with the help of
#'   \code{\link[randomcoloR]{distinctColorPalette}}, or this can be \code{NA}
#'   for no colors, or a vector of colors; the length of this vector
#'   of colors must match the number of bounded cells, which is displayed when
#'   you run the \code{\link{voronoi}} function and that you can also get by
#'   typing \code{attr(v, "nbounded")}
#' @param hue,luminosity if \code{colors = "random"}, these arguments are passed
#'   to \code{\link[randomcoloR]{randomColor}}
#' @param alpha opacity, a number between 0 and 1
#'   (used when \code{colors} is not \code{NA})
#' @param ... arguments passed to \code{\link{plotBoundedCell2D}} or
#'   \code{\link{plotBoundedCell3D}}
#'
#' @return No returned value.
#' @export
#' @importFrom randomcoloR randomColor distinctColorPalette
#' @importFrom scales alpha
#'
#' @note Sometimes, it is necessary to set the option \code{degenerate=TRUE}
#'   in the \code{\link{delaunay}} function in order to get a correct
#'   Voronoï diagram with the \code{plotVoronoiDiagram} function (I don't know
#'   why).
#'
#' @examples library(tessellation)
#' # 2D example: Fermat spiral
#' theta <- seq(0, 100, length.out = 300L)
#' x <- sqrt(theta) * cos(theta)
#' y <- sqrt(theta) * sin(theta)
#' pts <- cbind(x,y)
#' opar <- par(mar = c(0, 0, 0, 0), bg = "black")
#' # Here is a Fermat spiral:
#' plot(pts, asp = 1, xlab = NA, ylab = NA, axes = FALSE, pch = 19, col = "white")
#' # And here is its Voronoï diagram:
#' plot(NULL, asp = 1, xlim = c(-15, 15), ylim = c(-15, 15),
#'      xlab = NA, ylab = NA, axes = FALSE)
#' del <- delaunay(pts)
#' v <- voronoi(del)
#' length(Filter(isBoundedCell, v)) # 281 bounded cells
#' plotVoronoiDiagram(v, colors = viridisLite::turbo(281L))
#' par(opar)
#'
#' # 3D example: tetrahedron surrounded by three circles
#' tetrahedron <-
#'   rbind(
#'     c(2*sqrt(2)/3, 0, -1/3),
#'     c(-sqrt(2)/3, sqrt(2/3), -1/3),
#'     c(-sqrt(2)/3, -sqrt(2/3), -1/3),
#'     c(0, 0, 1)
#'   )
#' angles <- seq(0, 2*pi, length.out = 91)[-1]
#' R <- 2.5
#' circle1 <- t(vapply(angles, function(a) R*c(cos(a), sin(a), 0), numeric(3L)))
#' circle2 <- t(vapply(angles, function(a) R*c(cos(a), 0, sin(a)), numeric(3L)))
#' circle3 <- t(vapply(angles, function(a) R*c(0, cos(a), sin(a)), numeric(3L)))
#' circles <- rbind(circle1, circle2, circle3)
#' pts <- rbind(tetrahedron, circles)
#' d <- delaunay(pts, degenerate = TRUE)
#' v <- voronoi(d)
#' library(rgl)
#' open3d(windowRect = c(50, 50, 562, 562))
#' material3d(lwd = 2)
#' \donttest{plotVoronoiDiagram(v, luminosity = "bright")}
plotVoronoiDiagram <- function(
  v, colors = "random", hue = "random", luminosity = "light", alpha = 1, ...
){
  if(!inherits(v, "voronoi")){
    stop(
      "The argument `v` must be an output of the `voronoi` function.",
      call. = TRUE
    )
  }
  dimension <- length(v[[1L]][["site"]])
  cells <- Filter(isBoundedCell, v)
  ncells <- length(cells)
  if(ncells == 0L){
    stop(
      "This Vorono\u00ef tessellation has no bounded cells.",
      call. = TRUE
    )
  }
  if(identical(colors, "random")){
    colors <- scales::alpha(
      randomColor(ncells, hue = hue, luminosity = luminosity), alpha
    )
  }else if(identical(colors, "distinct")){
    colors <- scales::alpha(distinctColorPalette(ncells))
  }else if(identical(colors, NA)){
    colors <- rep(NA, ncells)
  }else{
    if(length(colors) != ncells){
      stop(
        sprintf(
          "There are %d bounded Vorono\u00ef cells and you supplied %d colors.",
          ncells, length(colors)
        )
      )
    }else{
      colors <- scales::alpha(colors)
    }
  }
  if(dimension == 2L){
    for(i in 1L:ncells){
      plotBoundedCell2D(
        cells[[i]], color = colors[i], check.bounded = FALSE, ...
      )
    }
  }else{
    for(i in 1L:ncells){
      plotBoundedCell3D(
        cells[[i]], facetsColor = colors[i], check.bounded = FALSE, ...
      )
    }
  }
}
