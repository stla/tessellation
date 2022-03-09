library(tessellation)
library(rgl)
library(hash)
library(Rvcg)

points <- rbind(
  c(0.5,0.5,0.5),
  c(0,0,0),
  c(0,0,1),
  c(0,1,0),
  c(0,1,1),
  c(1,0,0),
  c(1,0,1),
  c(1,1,0),
  c(1,1,1)
)
del <- delaunay(points, exteriorEdges = TRUE)

triangles <- vapply(del$tilefacets, function(x){
  as.integer(keys(x$subsimplex$vertices))
}, integer(3L))

mesh <- tmesh3d(
  t(points),
  triangles
)

shade3d(mesh, color = "orange", alpha=0.2)

tessellation:::exteriorDelaunayEdges(del)

