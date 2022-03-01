library(tessellation)

set.seed(666)
points <- matrix(rgamma(20L, 10, 1), nrow = 10L, ncol = 2L)

del <- delaunay(points)

vertices <- do.call(rbind, lapply(del[["vertices"]], `[[`, "point"))
triangles <- do.call(rbind, lapply(del[["tiles"]], `[[`, "facets"))
