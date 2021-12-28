cube <-
  rbind(
    c(-1, -1, -1),
    c( 1, -1, -1),
    c(-1,  1, -1),
    c( 1,  1, -1),
    c(-1, -1,  1),
    c( 1, -1,  1),
    c(-1,  1,  1),
    c( 1,  1,  1)
  )

xi_ <- seq(0, 2*pi, length.out = 91)[-1]
R <- 3
circle1 <- t(vapply(xi_, function(xi) R*c(cos(xi), sin(xi), 0), numeric(3L)))
circle2 <- t(vapply(xi_, function(xi) R*c(cos(xi), 0, sin(xi)), numeric(3L)))
circle3 <- t(vapply(xi_, function(xi) R*c(0, cos(xi), sin(xi)), numeric(3L)))
cube <- rbind(cube, circle1, circle2, circle3)

d <- delaunay(cube, degenerate = TRUE)
v <- voronoi(d)

library(paletteer)
library(rgl)
open3d(windowRect = c(50, 50, 512, 512))
bg3d("palegoldenrod")
plotVoronoiDiagram(v, colors = paletteer_c("grDevices::Dark 3", 8))

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
  frames = "zzzpic",
  convert = "echo \"%d %s %s %s\"",
  clean = FALSE
)

pngs <- list.files("./inst/ztemp/", pattern = "^zzzpic", full.names = TRUE)
library(gifski)
gifski(pngs, "CubeWithCircles.gif",
       width = 512, height = 512, delay = 1/10)
