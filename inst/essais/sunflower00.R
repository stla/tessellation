phi <- pi*(3 - sqrt(5))

n <- 1:150
c <- 0.1
r <- c * sqrt(n)
theta <- n * phi

x <- r * cos(theta)
y <- r * sin(theta)

opar <- par(mar = c(0, 0, 0, 0))
plot(x, y, pch = 19, asp = 1, axes = FALSE, xlab = NA, ylab = NA, cex = 0.5)

pts <- cbind(x, y)

d <- delaunay(pts)
plotDelaunay2D(d, asp = 1)
