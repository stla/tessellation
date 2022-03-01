library(tessellation)

f <- function(x, y){
  exp(-(x*x + y*y))
}
x <- seq(-5, 5, length.out = 100)
y <- seq(-5, 5, length.out = 100)
points <- as.matrix(transform( # data (x_i, y_i, z_i)
  expand.grid(x = x, y = y), z = f(x, y)
))

del <- delaunay(points, elevation = TRUE)

# an elevated Delaunay tessellation ####
f <- function(x, y){
  dnorm(x) * dnorm(y)
}
x <- y <- seq(-5, 5, length.out = 50)
grd <- expand.grid(x = x, y = y) # grid
points <- as.matrix(transform( # data (x_i, y_i, z_i)
  grd, z = f(x, y)
))
del <- delaunay(points, elevation = TRUE)
del[["volume"]] # close to 1, as expected
# plotting
library(rgl)
mesh <- del[["mesh"]]
open3d(windowRect = c(100, 100, 612, 356), zoom = 0.6)
aspect3d(1, 1, 20)
shade3d(mesh, color = "limegreen")
wire3d(mesh)
