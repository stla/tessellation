library(tessellation)
library(uniformly)
library(randomcoloR)

cube <- t(rgl::cube3d()$vb[-4L, ])

ptson <- runif_on_cube(20, 3)

pts <- rbind(cube, ptson)

d <- delaunay(pts)

sandwichedFacet <- function(tilefacet){
  length(tilefacet[["facetOf"]]) == 2L
}

x <- Filter(Negate(sandwichedFacet), d[["tilefacets"]])

sum(sapply(x, function(f) f$subsimplex$volume))
