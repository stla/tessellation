pts <- rbind(
  c(-5, -5,  16),
  c(-5,  8,   3),
  c(4,  -1,   3),
  c(4,  -5,   7),
  c(4,  -1, -10),
  c(4,  -5, -10),
  c(-5,  8, -10),
  c(-5, -5, -10)
)

dd <- delaunay(pts, atinfinity = FALSE, degenerate = FALSE)

vertices <- t(vapply(dd$vertices, `[[`, numeric(3), "point"))
simplicies <-
  lapply(lapply(dd[["tiles"]], `[[`, "simplex"), `[[`, "vertices")
edges <- unique(do.call(rbind, lapply(simplicies, function(simplex){
  t(combn(as.integer(keys(simplex)), 2L))
})))
simplex <- t(values(simplicies[[1]]))

triangles <- combn(4, 3)

library(rgl)
library(randomcoloR)
nsimplicies <- length(simplicies)
colors <- randomColor(nsimplicies, hue = "random", luminosity = "light")

for(i in seq_along(simplicies)){
  simplex <- t(values(simplicies[[i]]))
  for(j in 1:4){
    triangles3d(simplex[triangles[, j], ], color = colors[i], alpha = 0.3)
  }
}
for(i in 1:nrow(edges)){
  edge <- edges[i, ]
  p1 <- vertices[edge[1], ]
  p2 <- vertices[edge[2], ]
  lines3d(rbind(p1,p2), color = "black")
}
