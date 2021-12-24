library(tessellation)

# vertices ####
phi <- (1+sqrt(5))/2
a <- 1/sqrt(3)
b <- a/phi
c <- a*phi
vertices <-
  rbind(
    c( 0,  0,  0),
    c( a,  a,  a),
    c( a,  a, -a),
    c( a, -a,  a),
    c(-a, -a,  a),
    c(-a,  a, -a),
    c(-a,  a,  a),
    c( 0,  b, -c),
    c( 0, -b, -c),
    c( 0, -b,  c),
    c( c,  0, -b),
    c(-c,  0, -b),
    c(-c,  0,  b),
    c( b,  c,  0),
    c( b, -c,  0),
    c(-b, -c,  0),
    c(-b,  c,  0),
    c( 0,  b,  c),
    c( a, -a, -a),
    c( c,  0,  b),
    c(-a, -a, -a)
  )

d <- delaunay(vertices)

plotDelaunay3D(d, exteriorEdgesAsTubes = TRUE, tubeRadius = 0.03, tubeColor = "navy")

v <- voronoi(d)
open3d()
material3d(lwd = 2)
plotBoundedCell3D(v[[1]], facetsColor = "orange")
