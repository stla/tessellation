library(tessellation)
library(rgl)

approxUnitSphere <- subdivision3d(icosahedron3d(), depth = 3L)
approxUnitSphere$vb[4L, ] <-
  apply(approxUnitSphere$vb[1L:3L, ], 2L, function(x) sqrt(sum(x * x)))
homogeneousVertices <- t(approxUnitSphere$vb)
vertices <- asEuclidean(homogeneousVertices)

d <- delaunay(vertices, atinfinity = FALSE, degenerate = TRUE)

#volume(d); 4*pi/3

sandwichedFacet <- function(tilefacet){
  length(tilefacet[["facetOf"]]) == 2L
}

x <- Filter(Negate(sandwichedFacet), d[["tilefacets"]])

sum(sapply(x, function(f) f$subsimplex$volume))
