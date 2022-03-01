library(tessellation)
library(rgl)

cuboctahedron <- t(cuboctahedron3d()$vb[-4L, ])
xi_ <- seq(0, 2*pi, length.out = 91)[-1]
R <- 1.5
circle1 <- t(vapply(xi_, function(xi) R*c(cos(xi), sin(xi), 0), numeric(3L)))
circle2 <- t(vapply(xi_, function(xi) R*c(cos(xi), 0, sin(xi)), numeric(3L)))
circle3 <- t(vapply(xi_, function(xi) R*c(0, cos(xi), sin(xi)), numeric(3L)))
enclosedCuboctahedron <- rbind(cuboctahedron, circle1, circle2, circle3)

open3d(windowRect = c(50, 50, 562, 562))
shade3d(cuboctahedron3d(), color = "darkorange")
spheres3d(rbind(circle1, circle2, circle3), radius = 0.02)


d <- delaunay(enclosedCuboctahedron, degenerate = TRUE)
v <- voronoi(d)
#> Voronoï diagram with twelve bounded cells.


# make GIF ####
library(paletteer) # provides many color palettes
open3d(windowRect = c(50, 50, 562, 562))
bg3d("palegoldenrod")
plotVoronoiDiagram(v, colors = paletteer_c("grDevices::Dark 3", 12L))

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
  fps = 120,
  duration = 1,
  dir = ".",
  frames = "zzzpic",
  convert = "echo \"%d %s %s %s\"",
  clean = FALSE,
  webshot = FALSE
)

pngs <- list.files(".", pattern = "^zzzpic", full.names = TRUE)
library(gifski)
gifski(pngs, "voronoi_enclosed_cuboctahedron.gif",
       width = 512, height = 512, delay = 1/10)
file.remove(pngs)