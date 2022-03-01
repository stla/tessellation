library(tessellation)
library(rgl)
library(viridisLite)

xy <- as.matrix(expand.grid(x = 1L:nrow(volcano), y = 1L:ncol(volcano)))
z <-  c(volcano)
pts <- cbind(xy, z)

del <- delaunay(pts, elevation = TRUE)
mesh <- del[["mesh"]]

squaredNorms <- apply(mesh[["vb"]][-4L, ], 2L, crossprod)
normalizedSquaredNorms <- squaredNorms /  max(squaredNorms)

palette <- function(x) {
  RGB <- colorRamp(turbo(256))(x)
  rgb(RGB, maxColorValue = 255)
}

mesh$material <- list(color = palette(normalizedSquaredNorms))

open3d(windowRect = c(50, 50, 562, 306))
view3d(0, -40, zoom = 0.6)
aspect3d(5, 5, 1)
shade3d(mesh)
wire3d(mesh, color = "black")
