sameFamily <- function(tile1, tile2){
  family1 <- tile1[["family"]]
  family2 <- tile2[["family"]]
  if(is.na(family1) || is.na(family2)){
    return(FALSE)
  }
  family1 == family2
}

edgesFromTileFacet <- function(tessellation, tilefacet){
  tileindices <- tilefacet[["facetOf"]]
  tiles <- tessellation[["tiles"]]
  tile1 <- tiles[[tileindices[1L]]]
  c1 <- tile1[["simplex"]][["circumcenter"]]
  if(length(tileindices) == 1L){
    return(IEdge$new(c1, tilefacet[["normal"]]))
  }
  tile2 <- tiles[[tileindices[2L]]]
  c2 <- tile2[["simplex"]][["circumcenter"]]
  if(sameFamily(tile1, tile2) || isTRUE(all.equal(c1, c2))){
    return(NULL)
  }
  Edge$new(c1, c2)
}

voronoiCell <- function(facetsQuotienter, edgeTransformer){
  function(tessellation, vertexId){
    tilefacets <- facetsQuotienter(
      unname(vertexNeighborFacets(tessellation, vertexId))
    )
    edges <- Filter(Negate(is.null), lapply(tilefacets, function(tilefacet){
      edgesFromTileFacet(tessellation, tilefacet)
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
  L <- lapply(1:nrow(vertices), function(i) cellGetter(tessellation, i))
  zip(vertices, L)
}

#' @title Voronoï tessellation
#' @description Voronoï tessellation from Delaunay tessellation; this is a list
#'   of pairs made of a site (a vertex) and a list of edges.
#'
#' @param tessellation output of \code{\link{delaunay}}
#'
#' @return A list representing the Voronoï tessellation.
#' @export
#'
#' @examples library(tessellation)
#' d <- delaunay(centricCuboctahedron())
#' voronoi(d)
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
  all(vapply(cell, function(edge) inherits(edge, "Edge"), logical(1L)))
}

#' @title Vertices of a bounded cell
#' @description Get all vertices of a bounded cell, without duplicates.
#'
#' @param cell A bounded Voronoï cell.
#'
#' @return A matrix, each row represents a vertex.
#' @export
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
