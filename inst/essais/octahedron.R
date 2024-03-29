library(tessellation)
library(rgl)

octahedron <- t(octahedron3d()$vb[-4, ])
square1 <- rbind(
  c(3, 3, 0),
  c(3, -3, 0),
  c(-3, 3, 0),
  c(-3, -3, 0)
)
square2 <- rbind(
  c(4, 0, 4),
  c(4, 0, -4),
  c(-4, 0, 4),
  c(-4, 0, -4)
)
square3 <- rbind(
  c(0, 5, 5),
  c(0, 5, -5),
  c(0, -5, 5),
  c(0, -5, -5)
)

pts <- rbind(octahedron, square1, square2, square3)

del <- delaunay(pts)

v <- voronoi(del)

plotVoronoiDiagram(v, colors = rainbow(6))

