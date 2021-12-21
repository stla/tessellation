library(tessellation)
library(uniformly)
library(randomcoloR)

tetrahedron <- rbind(
  c(sqrt(3)/3, 0 ,0),
  c(-sqrt(3)/6, 0.5, 0),
  c(-sqrt(3)/6, -0.5, 0),
  c(0, 0, sqrt(6)/3)
)

ptsin <- runif_in_simplex(10, tetrahedron)

pts <- rbind(tetrahedron, ptsin)

d <- delaunay(pts)

v <- voronoi(d)

ncells <- length(v)

colors <- distinctColorPalette(ncells)

open3d(windowRect = c(50, 50, 562, 562))
for(i in 1:ncells){
  if(isBoundedCell(v[[i]]))
    plotBoundedCell(v[[i]][["cell"]], facetsColor = colors[i])
}

