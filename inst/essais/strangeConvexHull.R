library(tessellation)

alpha <- seq(0, 2*pi, length.out = 150)[-1]
points <- t(vapply(alpha, function(x){
  c(
    sin(x)*cos(2*x), sin(x)*sin(2*x), cos(x)
  )
}, numeric(3L)))
points <- rbind(c(0, 0, 0), points)


d <- delaunay(points, degenerate = TRUE)
v <- voronoi(d)

library(rgl)
open3d(windowRect = c(50, 50, 450, 450))
plotBoundedCell3D(v[[1]])


# animation ####
M <- par3d("userMatrix")
movie3d(
  par3dinterp(
    time = seq(0, 1, len = 9),
    userMatrix = list(
      M,
      rotate3d(M, pi, 1, 0, 0),
      rotate3d(M, pi, 1, 1, 0),
      rotate3d(M, pi, 1, 1, 1),
      rotate3d(M, pi, 0, 1, 1),
      rotate3d(M, pi, 0, 1, 0),
      rotate3d(M, pi, 1, 0, 1),
      rotate3d(M, pi, 0, 0, 1),
      M
    )
  ),
  fps = 100,
  duration = 1,
  dir = "./inst/ztemp/",
  frames = "oooooo",
  convert = "echo \"%d %s %s %s\"",
  clean = FALSE
)

pngs <- list.files("./inst/ztemp/", pattern = "^oooooo", full.names = TRUE)
library(gifski)
gifski(pngs, "strangeVoronoiCell_400x.gif",
       width = 400, height = 400, delay = 1/10)
