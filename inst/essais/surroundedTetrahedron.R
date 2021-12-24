library(RColorBrewer)
library(tessellation)
library(rgl)

tetrahedron <-
  rbind(
    c(2*sqrt(2)/3, 0, -1/3),
    c(-sqrt(2)/3, sqrt(2/3), -1/3),
    c(-sqrt(2)/3, -sqrt(2/3), -1/3),
    c(0, 0, 1)
  )

faces <- combn(4L, 3L)
edges <- combn(4L, 2L)

open3d(windowRect = c(50, 50, 562, 562))
for(j in 1L:ncol(faces)){
  triangles3d(tetrahedron[faces[, j], ], color = "darkred")
}
for(j in 1L:ncol(edges)){
  shade3d(
    cylinder3d(tetrahedron[edges[, j], ], sides = 60, radius = 0.03),
    color = "yellow"
  )
}
spheres3d(tetrahedron, radius = 0.05, color = "yellow")

alpha <- seq(0, 2*pi, length.out = 91)[-1]
R <- 2.5
circle1 <- t(vapply(alpha, function(al) R*c(cos(al), sin(al), 0), numeric(3L)))
circle2 <- t(vapply(alpha, function(al) R*c(cos(al), 0, sin(al)), numeric(3L)))
circle3 <- t(vapply(alpha, function(al) R*c(0, cos(al), sin(al)), numeric(3L)))
circles <- rbind(circle1, circle2, circle3)
spheres3d(circles, radius = 0.04, color = "black")

object <- rbind(tetrahedron, circles)

colors <- brewer.pal(n = 4, name = "Dark2")
d <- delaunay(object, degenerate = TRUE)
v <- voronoi(d)
library(rgl)
open3d(windowRect = c(50, 50, 562, 562))
for(i in 1:4){
  plotBoundedCell(v[[i]], facetsColor = colors[i])
}

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
  frames = "zpictetra",
  convert = "echo \"%d %s %s %s\"",
  clean = FALSE
)

pngs <- list.files("./inst/ztemp/", pattern = "^zpictetra", full.names = TRUE)
gifski(pngs, "surroundedTetrahedron.gif",
       width = 512, height = 512, delay = 1/12)



###########


tetrahedron <-
  rbind(
    c(2*sqrt(2)/3, 0, -1/3),
    c(-sqrt(2)/3, sqrt(2/3), -1/3),
    c(-sqrt(2)/3, -sqrt(2/3), -1/3),
    c(0, 0, 1)
  )
angles <- seq(0, 2*pi, length.out = 91)[-1]
R <- 2.5
circle1 <- t(vapply(angles, function(a) R*c(cos(a), sin(a), 0), numeric(3L)))
circle2 <- t(vapply(angles, function(a) R*c(cos(a), 0, sin(a)), numeric(3L)))
circle3 <- t(vapply(angles, function(a) R*c(0, cos(a), sin(a)), numeric(3L)))
circles <- rbind(circle1, circle2, circle3)
pts <- rbind(tetrahedron, circles)
d <- delaunay(pts, degenerate = TRUE)
v <- voronoi(d)
library(rgl)
open3d(windowRect = c(50, 50, 562, 562))
plotVoronoiDiagram(v, luminosity = "dark")
