library(viridisLite)
library(randomcoloR)
library(tessellation)
library(rgl)

dodecahedron <-
  t(dodecahedron3d()$vb[-4,])

alpha <- seq(0, 2*pi, length.out = 91)[-1]
R <- 3
circle1 <- t(vapply(alpha, function(al) R*c(cos(al), sin(al), 0), numeric(3L)))
circle2 <- t(vapply(alpha, function(al) R*c(cos(al), 0, sin(al)), numeric(3L)))
circle3 <- t(vapply(alpha, function(al) R*c(0, cos(al), sin(al)), numeric(3L)))
dodecahedron <- rbind(dodecahedron, circle1, circle2, circle3)


colors <- distinctColorPalette(20)
d <- delaunay(dodecahedron, degenerate = FALSE)
v <- voronoi(d)
library(rgl)
open3d(windowRect = c(50, 50, 562, 562))
for(i in 1:20){
  if(isBoundedCell(v[[i]][["cell"]]))
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
  frames = "zzzpicdodec",
  convert = "echo \"%d %s %s %s\"",
  clean = FALSE
)
