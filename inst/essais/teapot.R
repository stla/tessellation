utah <- rbind(teapot(), c(0, 1.5, 0))

d <- delaunay(utah)
v <- voronoi(d)
cell1977 <- v[[1977]]
cellvertices <- cellVertices(cell1977)
# computes the triangulated convex hull:
library(cxhull)
h <- cxhull(cellvertices, triangulate = TRUE)
# plot:
open3d(windowRect = c(50,50,450,450))
material3d(specular = "blue")
view3d(10, 80, zoom = 0.75)
for(i in 1:length(h$facets)){
  triangle <- t(sapply(h$facets[[i]]$vertices,
                       function(id) h$vertices[[as.character(id)]]$point))
  triangles3d(triangle, color = "blue")
}

invisible(lapply(cell1977[["cell"]], function(edge) edge$plot()))
