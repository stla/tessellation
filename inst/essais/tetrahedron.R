#library(tessellation)

tetrahedron <-
  rbind(
    c(2*sqrt(2)/3, 0, -1/3),
    c(-sqrt(2)/3, sqrt(2/3), -1/3),
    c(-sqrt(2)/3, -sqrt(2/3), -1/3),
    c(0, 0, 1)
  )

library(uniformly)
set.seed(314)
randomPoints <- runif_in_tetrahedron(
  3, tetrahedron[1, ], tetrahedron[2, ], tetrahedron[3, ], tetrahedron[4, ]
)

library(rgl)
plotTetrahedron <- function(tetrahedron, alpha){
  faces <- combn(4L, 3L)
  for(j in 1L:ncol(faces)){
    triangles3d(tetrahedron[faces[, j], ], color = "orange", alpha = alpha)
  }
  edges <- combn(4L, 2L)
  for(j in 1L:ncol(edges)){
    shade3d(
      cylinder3d(tetrahedron[edges[, j], ], sides = 60, radius = 0.03),
      color = "yellow"
    )
  }
  spheres3d(tetrahedron, radius = 0.05, color = "yellow")
}

open3d(windowRect = c(50, 50, 562, 562))
bg3d("slategray")
plotTetrahedron(tetrahedron, alpha = 0.3)
spheres3d(randomPoints, radius = 0.04, color = "black")

# movie3d(
#   spin3d(axis = c(0, 0, 1), rpm = 12),
#   duration = 5, fps = 12,
#   movie = "Tetrahedron", dir = ".",
#   convert = "magick convert -dispose previous -loop 0 -delay 1x%d %s*.png %s.%s",
#   startTime = 1/60
# )

library(tessellation)
pts <- rbind(tetrahedron, randomPoints)
del <- delaunay(pts)
open3d(windowRect = c(50, 50, 562, 562))
bg3d("slategray")
material3d(lwd = 2)
plotDelaunay3D(del, color=TRUE, luminosity = "bright", alpha=0.5)

# movie3d(
#   spin3d(axis = c(0, 0, 1), rpm = 12),
#   duration = 5, fps = 12,
#   movie = "TetrahedronDelaunay", dir = ".",
#   convert = "magick convert -dispose previous -loop 0 -delay 1x%d %s*.png %s.%s",
#   startTime = 1/60
# )




tiles <- getDelaunaySimplicies(del)
for(i in seq_along(tiles)){
  open3d(windowRect = c(50, 50, 306, 306), zoom = 0.8)
  bg3d("slategray")
  material3d(lwd = 2)
  plotDelaunay3D(del, color=FALSE)
  plotTetrahedron(tiles[[i]], alpha = 0.6)
  rgl.snapshot(sprintf("TetrahedronTile%02d.png", i))
  close3d()
}



del <- delaunay(pts, exteriorEdges = TRUE)
open3d(windowRect = c(50, 50, 562, 562))
bg3d("slategray")
material3d(lwd = 2)
plotDelaunay3D(
  del, color=TRUE, luminosity = "bright", alpha = 0.2,
  exteriorEdgesAsTubes = TRUE, tubeRadius = 0.03, tubeColor = "navy"
)
# movie3d(
#   spin3d(axis = c(0, 0, 1), rpm = 12),
#   duration = 5, fps = 12,
#   movie = "TetrahedronDelaunayWithTubes", dir = ".",
#   convert = "magick convert -dispose previous -loop 0 -delay 1x%d %s*.png %s.%s",
#   startTime = 1/60
# )

# Voronoï ####
v <- voronoi(del)
open3d(windowRect = c(50, 50, 562, 562))
bg3d("palegoldenrod")
material3d(lwd = 2)
plotVoronoiDiagram(v, luminosity = "dark")

# # animation ####
# M <- par3d("userMatrix")
# movie3d(
#   par3dinterp(
#     time = seq(0, 1, len = 9),
#     userMatrix = list(
#       M,
#       rotate3d(M, pi, 1, 0, 0),
#       rotate3d(M, pi, 1, 1, 0),
#       rotate3d(M, pi, 1, 1, 1),
#       rotate3d(M, pi, 0, 1, 1),
#       rotate3d(M, pi, 0, 1, 0),
#       rotate3d(M, pi, 1, 0, 1),
#       rotate3d(M, pi, 0, 0, 1),
#       M
#     )
#   ),
#   fps = 100,
#   duration = 1,
#   dir = ".",
#   frames = "zVoronoiTetrahedron",
#   convert = "echo \"%d %s %s %s\"",
#   clean = FALSE
# )
#
# pngs <- list.files(".", pattern = "^zVoronoiTetrahedron", full.names = TRUE)
# library(gifski)
# gifski(pngs, "VoronoiTetrahedron.gif",
#        width = 512, height = 512, delay = 1/10)
#####################

xi_ <- seq(0, 2*pi, length.out = 91)[-1]
R <- 2
circle1 <- t(vapply(xi_, function(xi) R * c(cos(xi), sin(xi), 0), numeric(3L)))
circle2 <- t(vapply(xi_, function(xi) R * c(cos(xi), 0, sin(xi)), numeric(3L)))
circles <- rbind(circle1, circle2)

open3d(windowRect = c(50, 50, 562, 562), zoom=0.7)
bg3d("palegoldenrod")
plotTetrahedron(tetrahedron, alpha = 0.3)
spheres3d(randomPoints, radius = 0.04, color = "black")
spheres3d(circles, radius = 0.04, color = "black")


object <- rbind(tetrahedron, randomPoints, circles)
del <- delaunay(object)
v <- voronoi(del)
library(viridisLite)
open3d(windowRect = c(50, 50, 562, 562))
bg3d("palegoldenrod")
plotVoronoiDiagram(v, colors = turbo(7))
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
  frames = "zVoronoi",
  convert = "echo \"%d %s %s %s\"",
  clean = FALSE
)

pngs <- list.files(".", pattern = "^zVoronoi", full.names = TRUE)
library(gifski)
gifski(pngs, "VoronoiTetrahedronWithCircles.gif",
       width = 512, height = 512, delay = 1/10)

