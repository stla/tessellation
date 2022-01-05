library(tessellation)
library(rgl)

approxUnitSphere <- subdivision3d(icosahedron3d(), depth = 5L)
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


library(misc3d)
fx <- function(theta, phi){
  cos(theta) * sin(phi)
}
fy <- function(theta, phi){
  sin(theta) * sin(phi)
}
fz <- function(theta, phi){
  cos(phi)
}

sphere <- parametric3d(fx, fy, fz, umin = 0, umax = 2*pi, vmin = 0, vmax = pi,
                       engine = "rgl", n = 30, color = "black", fill = FALSE)

sphere <- parametric3d(fx, fy, fz, umin = 0, umax = 2*pi, vmin = 0, vmax = pi,
                       engine = "none", n = 40)

vertices <- rbind(sphere$v1, sphere$v2, sphere$v3)

d <- delaunay(vertices)#, degenerate = TRUE)
volume(d)  # 4*pi/3
surface(d) # 4*pi

