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
    return(IEdge$new(c1, tilefacet[["subsimplex"]][["normal"]]))
  }
  tile2 <- tiles[[tileindices[2L]]]
  c2 <- tile2[["simplex"]][["circumcenter"]]
  if(sameFamily(tile1, tile2) || isTRUE(all.equal(c1, c2))){
    return(NULL)
  }
  Edge$new(center1, center2)
}

voronoiCell <- function(facetsQuotienter, edgeTransformer, tessellation, vertexId){
  tilefacets <- facetsQuotienter(vertexNeighborFacets(tessellation, vertexId))
  edges <- Filter(Negate(is.null), lapply(tilefacets, function(tilefacet){
    edgesFromTileFacet(tessellation, tilefacet)
  }))
  edgeTransformer(edges)
}

