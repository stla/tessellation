utah <- rbind(teapot(), c(0, 1.5, 0))

d <- delaunay(utah)
v <- voronoi(d)
cell1977 <- v[[1977]]

plotBoundedCell3D(cell1977, facetsColor = "navy")

pts <- cellVertices(cell1977)
pts <- rbind(pts, c(0, 1.5, 0))
d <- delaunay(pts, degenerate = FALSE)
v <- voronoi(d)
plotBoundedCell3D(v[[803]], facetsColor = "navy")

plotVoronoiDiagram(v)
