phi <- pi*(3 - sqrt(5))

n <- c(0, 50:150)
c <- 1
r <- c * n*n
theta <- n * phi


x <- r * cos(theta)
y <- r * sin(theta)

opar <- par(mar = c(0, 0, 0, 0))
plot(x, y, pch = 19, asp = 1, axes = FALSE, xlab = NA, ylab = NA, cex = 0.5)

pts <- cbind(x, y)

d <- delaunay(pts)
plotDelaunay2D(d, asp = 1)



angles <- seq(0, 2*pi, length.out=91)[-1]
circle <- 25000 * cbind(cos(angles), sin(angles))
pts <- rbind(pts, circle)
d <- delaunay(pts)
v <- voronoi(d)
opar <- par(mar = c(0, 0, 0, 0))
plot(pts, type = "n", xlab = NA, ylab = NA, asp = 1, axes = FALSE)
plotVoronoiDiagram(
  v, luminosity = "dark"
)
par(opar)
