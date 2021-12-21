cube <-
  rbind(
    c(-1, -1, -1),
    c(1, -1, -1),
    c(-1, 1, -1),
    c(1, 1, -1),
    c(-1, -1, 1),
    c(1, -1, 1),
    c(-1, 1, 1),
    c(1, 1, 1)
  )
alpha <- seq(0, 2*pi, length.out = 91)[-1]
R <- 3
circle1 <- t(vapply(alpha, function(al) R*c(cos(al), sin(al), 0), numeric(3L)))
circle2 <- t(vapply(alpha, function(al) R*c(cos(al), 0, sin(al)), numeric(3L)))
circle3 <- t(vapply(alpha, function(al) R*c(0, cos(al), sin(al)), numeric(3L)))
cube <- rbind(cube, circle1, circle2, circle3)


d <- delaunay(cube, degenerate = TRUE)
v <- voronoi(d)
library(rgl)
open3d(windowRect = c(50, 50, 562, 562))
for(i in 1:8){
  #cell000 <- Filter(function(x) inherits(x, "Edge"), v[[i]][["cell"]])
  lapply(v[[i]][["cell"]], function(edge) edge$plot())
}
