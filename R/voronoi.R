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
  if(length(tileindices) == 1L){
    return(newIEdge(c1, tilefacet[["normal"]]))
  }
  tile2 <- tiles[[tileindices[2L]]]
  c2 <- tile2[["simplex"]][["circumcenter"]]
  if(sameFamily(tile1, tile2) || isTRUE(all.equal(c1, c2))){
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
    edgeTransformer(edges)
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
  L <- lapply(1L:nrow(vertices), function(i) cellGetter(tessellation, i))
  zip(vertices, L)
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
#'   \code{\link{plotBoundedCell}}
#' @export
#'
#' @examples library(tessellation)
#' d <- delaunay(centricCuboctahedron())
#' v <- voronoi(d)
#' # the Voronoï diagram has 13 cells (one for each site):
#' length(v)
#' # there is only one bounded cell:
#' length(Filter(isBoundedCell, v))
voronoi <- function(tessellation){
  voronoi0(voronoiCell(identity, identity), tessellation)
}

#' @title Is this cell bounded?
#' @description Check whether a Voronoï cell is bounded, i.e. contains only
#'   finite edges.
#'
#' @param cell Voronoï cell
#'
#' @return A Boolean value, whether the cell is bounded.
#' @export
#' @importFrom sets is.tuple
isBoundedCell <- function(cell){
  if(is.tuple(cell)){
    cell <- cell[["cell"]]
  }
  all(vapply(cell, function(edge){
    inherits(edge, "Edge2") || inherits(edge, "Edge3")
  }, logical(1L)))
}

#' @title Vertices of a bounded cell
#' @description Get all vertices of a bounded cell, without duplicates.
#'
#' @param cell a bounded Voronoï cell
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
cellVertices <- function(cell){
  if(!isBoundedCell(cell)){
    stop(
      "This function applies to bounded cells only.",
      call. = TRUE
    )
  }
  if(is.tuple(cell)){
    cell <- cell[["cell"]]
  }
  stackedVertices <- lapply(cell, function(edge) edge$stack())
  unique(do.call(rbind, stackedVertices))
}

#' @title Plot a bounded Voronoï 3D cell
#' @description Plot a bounded Voronoï 3D cell with \strong{rgl}.
#'
#' @param cell a bounded Voronoï 3D cell
#' @param edgesAsTubes Boolean, whether to plot edges as tubes or as lines
#' @param tubeRadius radius of the tubes if \code{edgesAsTubes = TRUE}
#' @param tubeColor color of the tubes if \code{edgesAsTubes = TRUE}
#' @param facetsColor color of the facets; \code{NULL} for no color
#' @param alpha opacity of the facets, a number between 0 and 1
#'
#' @return No value, this function just plots the cell.
#' @export
#' @importFrom cxhull cxhull
#' @importFrom rgl lines3d cylinder3d shade3d
#' @importFrom sets is.tuple
#'
#' @examples library(tessellation)
#' d <- delaunay(centricCuboctahedron())
#' v <- voronoi(d)
#' cell13 <- v[[13]]
#' isBoundedCell(cell13) # TRUE
#' library(rgl)
#' open3d(windowRect = c(50, 50, 562, 562))
#' plotBoundedCell3D(
#'   cell13, edgesAsTubes = TRUE, tubeRadius = 0.03, tubeColor = "yellow",
#'   facetsColor = "navy", alpha = 0.7
#' )
plotBoundedCell3D <- function(
  cell, edgesAsTubes = FALSE, tubeRadius, tubeColor,
  facetsColor = NULL, alpha = 1
){
  if(!isBoundedCell(cell)){
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
  if(edgesAsTubes || !is.null(facetsColor)){
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
  if(!is.null(facetsColor)){
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
#' @param tubeColor color of the tubes if \code{edgesAsTubes = TRUE}
#' @param border color of the borders of the cell; \code{NA} for no color
#' @param color color of the cell; \code{NA} for no color
#' @param ... graphical parameters for the borders
#'
#' @return No value, this function just plots the cell.
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
#' plot(centricSquare, type = "n", asp = 1)
#' plotBoundedCell2D(cell5, color = "pink")
plotBoundedCell2D <- function(
  cell, border = "black", color = NA, ...
){
  if(!isBoundedCell(cell)){
    stop(
      "This function applies to bounded cells only.",
      call. = TRUE
    )
  }
  if(is.tuple(cell)){
    cell <- cell[["cell"]]
  }
  if(!is.na(color)){
    stackedVertices <- lapply(cell, function(edge) edge$stack())
    vertices <- unique(do.call(rbind, stackedVertices))
    polygon(vertices, border = NA, col = color)
  }
  if(!is.na(border)){
    for(edge in cell){
      edge$plot(color = border, ...)
    }
  }
}
