centricCuboctahedron <-
  rbind(
    c(-1, -1, 0),
    c(-1, 1, 0),
    c(1, -1, 0),
    c(1, 1, 0),
    c(-1, 0, -1),
    c(-1, 0, 1),
    c(1, 0, -1),
    c(1, 0, 1),
    c(0, -1, -1),
    c(0, -1, 1),
    c(0, 1, -1),
    c(0, 1, 1),
    c(0, 0, 0)
  )
d <- delaunay(centricCuboctahedron)
v <- voronoi(d)
cell000 <- v[[13]]
library(rgl)
open3d(windowRect = c(50, 50, 562, 562))
lapply(cell000[[2]], function(edge) edge$plot())
