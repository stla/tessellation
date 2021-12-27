
xi_ <- seq(0, 2*pi, length.out = 91)[-1]
R <- 3
circle1 <- t(vapply(xi_, function(xi) R * c(cos(xi), sin(xi), 0), numeric(3L)))
circle2 <- t(vapply(xi_, function(xi) R * c(cos(xi), 0, sin(xi)), numeric(3L)))
circle3 <- t(vapply(xi_, function(xi) R * c(0, cos(xi), sin(xi)), numeric(3L)))
circles <- rbind(circle1, circle2, circle3)

dodecahedron <- t(dodecahedron3d()$vb[-4,])

pts <- rbind(dodecahedron, circles)

del <- delaunay(pts, degenerate = TRUE)
v <- voronoi(del)

open3d(windowRect = c(50, 50, 562, 562))
bg3d("aliceblue")
plotVoronoiDiagram(v, colors = rainbow(20))

movie3d(
  spin3d(axis = c(0, 0, 1), rpm = 12),
  duration = 5, fps = 12,
  movie = "DodecahedronWithCircles", dir = ".",
  convert = "magick convert -dispose previous -loop 0 -delay 1x%d %s*.png %s.%s",
  startTime = 1/60
)







C0 <- (4 - sqrt(2)) / 7
C1 <- sqrt(2)
vertices <- rbind(
  c( 0.0,  0.0,   C1),
  c( 0.0,  0.0,  -C1),
  c(  C1,  0.0,  0.0),
  c( -C1,  0.0,  0.0),
  c( 0.0,   C1,  0.0),
  c( 0.0,  -C1,  0.0),
  c(-1.0,  0.0, -1.0),
  c(-1.0,  0.0,  1.0),
  c( 1.0,  0.0, -1.0),
  c( 1.0,  0.0,  1.0),
  c(-1.0, -1.0,  0.0),
  c(-1.0,  1.0,  0.0),
  c( 1.0, -1.0,  0.0),
  c( 1.0,  1.0,  0.0),
  c( 0.0, -1.0, -1.0),
  c( 0.0, -1.0,  1.0),
  c( 0.0,  1.0, -1.0),
  c( 0.0,  1.0,  1.0),
  c( -C0,  -C0,  -C0),
  c( -C0,  -C0,   C0),
  c( -C0,   C0,  -C0),
  c( -C0,   C0,   C0),
  c(  C0,  -C0,  -C0),
  c(  C0,  -C0,   C0),
  c(  C0,   C0,  -C0),
  c(  C0,   C0,   C0)
)

pts <- rbind(vertices, circles)


del <- delaunay(pts, degenerate = TRUE)
v <- voronoi(del)

open3d(windowRect = c(50, 50, 562, 562))
bg3d("aliceblue")
plotVoronoiDiagram(v, colors = plasma(26))

