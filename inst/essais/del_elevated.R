library(tessellation)

f <- function(x, y){
  exp(-(x*x + y*y))
}
x <- seq(0, 5, length.out = 30)
y <- seq(-5, 0, length.out = 30)
points <- as.matrix(transform( # data (x_i, y_i, z_i)
  expand.grid(x = x, y = y), z = f(x, y)
))

del <- delaunay(points, elevation = TRUE)
dd <- RCGAL::delaunay(points, elevation = TRUE)
del$volume
dd$volume

library(rgl)
mesh <- del$mesh
open3d(windowRect = c(50, 50, 562, 562))
shade3d(mesh, color = "yellow")
wire3d(mesh, color = "black")


xy <- as.matrix(expand.grid(x = x, y = y))
# xyo <- xy[order(rowSums(xy)), ]
# xy <- xy[order(xy[, 2L]-xy[, 1L]), ]
xy <- xy[order(round(rowSums(xy),6), xy[, 2L]-xy[, 1L]), ]

library(RTriangle)
tt <- triangulate(pslg(xy))

library(tripack)
plot(tri.mesh(xy[, 1L], xy[, 2L]), do.points = FALSE)

library(rgl)
open3d()
aspect3d(1,1,5)
wire3d(del$mesh, color="blue")
open3d()
aspect3d(1,1,5)
wire3d(dd$mesh, color="black")

pts2 <- matrix(NA_real_, nrow = 100, ncol = 2)
counter <- 1
for(i in 1:5){
  for(j in 1:5){
    k <- (i+j-2)*(i+j-1)/2+i
    print(k)
    pts2[counter, ] <- pts[k, ]
    counter <- counter+1
  }
}

pts2 <- matrix(NA_real_, nrow = 100, ncol = 2)
counter <- 1
for(j in 1:5){
  for(i in 1:5){
    k <- (i+j-2)*(i+j-1)/2+j
    print(k)
    #pts2[counter, ] <- pts[k, ]
    counter <- counter+1
  }
}

compositions(2, 2, FALSE)
compositions(3, 2, FALSE)
compositions(4, 2, FALSE)

f <- function(n){
  a <- partitions::compositions(n, 2,FALSE)
  a[, apply(a, 2, max) <= 20, drop=FALSE]
}
sum(sapply(2:20, function(n) ncol(f(n))))

pts <- matrix(NA_real_, nrow = 0, ncol = 2)
for(n in 2:40){
  a <- f(n)
#  if(n%%2 == 1) a <- a[, rev(seq_len(ncol(a)))]
  for(i in 1:ncol(a)){
    col <- a[, i]
    pts <- rbind(pts, c(x[col[1]], y[col[2]]))
  }
}

del <- tessellation::delaunay(pts)
tessellation::plotDelaunay2D(del)

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
