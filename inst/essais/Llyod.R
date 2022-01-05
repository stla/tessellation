library(tessellation)
library(uniformly)

square <- rbind(
  c(10, 10),
  c(-10, 10),
  c(-10, -10),
  c(10, -10)
)

pts <- NULL
for(i in 1:7){
  for(j in 1:7){
    pts <- rbind(pts, c(as.numeric(i), as.numeric(j)))
  }
}

set.seed(666)
pts <- runif_in_cube(50, 2, r = 100)
plot(pts, asp = 1, pch = 19)

toRemove <- NULL
D <- as.matrix(dist(pts))
for(i in 2:nrow(D)){
  row <- D[i, 1:(i-1)]
  if(any(row < 12)){
    toRemove <- c(toRemove, i)
  }
}
plot(pts[-toRemove,], asp = 1, pch = 19)



d <- delaunay(pts)
plotDelaunay2D(d)

tiles <- d$tiles
pts <- NULL
for(tile in tiles){
  print(tile$simplex$circumcenter)
  pts <- rbind(pts, rowMeans(values(tile$simplex$vertices)))
}

plot(pts, asp = 1, pch = 19)

v <- voronoi(d)
plotVoronoiDiagram(v)

