
rmesh <- Rvcg::vcgPlyRead("septuaginta.ply")

library(rgl) # provides `cuboctahedron3d()`
septuaginta_vertices <- t(rmesh$vb[-4L, ])
xi_ <- seq(0, 2*pi, length.out = 91)[-1L]
R <- 1.5
circle1 <- t(vapply(xi_, function(xi) R*c(cos(xi), sin(xi), 0), numeric(3L)))
circle2 <- t(vapply(xi_, function(xi) R*c(cos(xi), 0, sin(xi)), numeric(3L)))
circle3 <- t(vapply(xi_, function(xi) R*c(0, cos(xi), sin(xi)), numeric(3L)))
enclosedSeptuaginta <-
  rbind(septuaginta_vertices, circle1, circle2, circle3)

open3d(windowRect = c(50, 50, 562, 562))
view3d(20, zoom = 0.65)
shade3d(rmesh, color = "darkorange", shininess = 10)
wire3d(rmesh, color = "darkslategray4", lwd = 4)
spheres3d(rbind(circle1, circle2, circle3), radius = 0.04)

library(tessellation)
del <- delaunay(enclosedSeptuaginta, degenerate = TRUE)
v <- voronoi(del)

library(paletteer) # provides many color palettes
open3d(windowRect = c(50, 50, 562, 562))
bg3d("palegoldenrod")
plotVoronoiDiagram(v, colors = paletteer_c("grDevices::Dark 3", 62L))


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
  dir = ".",
  frames = "zzpic",
  convert = FALSE,
  webshot = FALSE,
  clean = FALSE
)

library(gifski)
pngs <- Sys.glob("zzpic*.png")
gifski(
  pngs,
  "septuaginta.gif",
  width = 512, height = 512,
  delay = 1/8
)
