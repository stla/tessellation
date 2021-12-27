library(tessellation)

sqrt2plus1 = sqrt(2) + 1

# vertices ####################################################################
vertices = rbind(
  c(-1, -1, -1, -sqrt2plus1),
  c(1, -1, -1, -sqrt2plus1),
  c(-1, 1, -1, -sqrt2plus1),
  c(1, 1, -1, -sqrt2plus1),
  c(-1, -1, 1, -sqrt2plus1),
  c(1, -1, 1, -sqrt2plus1),
  c(-1, 1, 1, -sqrt2plus1),
  c(1, 1, 1, -sqrt2plus1),
  c(-1, -1, -1, sqrt2plus1),
  c(1, -1, -1, sqrt2plus1),
  c(-1, 1, -1, sqrt2plus1),
  c(1, 1, -1, sqrt2plus1),
  c(-1, -1, 1, sqrt2plus1),
  c(1, -1, 1, sqrt2plus1),
  c(-1, 1, 1, sqrt2plus1),
  c(1, 1, 1, sqrt2plus1),
  c(-1, -1, -sqrt2plus1, -1),
  c(1, -1, -sqrt2plus1, -1),
  c(-1, 1, -sqrt2plus1, -1),
  c(1, 1, -sqrt2plus1, -1),
  c(-1, -1, sqrt2plus1, -1),
  c(1, -1, sqrt2plus1, -1),
  c(-1, 1, sqrt2plus1, -1),
  c(1, 1, sqrt2plus1, -1),
  c(-1, -1, -sqrt2plus1, 1),
  c(1, -1, -sqrt2plus1, 1),
  c(-1, 1, -sqrt2plus1, 1),
  c(1, 1, -sqrt2plus1, 1),
  c(-1, -1, sqrt2plus1, 1),
  c(1, -1, sqrt2plus1, 1),
  c(-1, 1, sqrt2plus1, 1),
  c(1, 1, sqrt2plus1, 1),
  c(-1, -sqrt2plus1, -1, -1),
  c(1, -sqrt2plus1, -1, -1),
  c(-1, sqrt2plus1, -1, -1),
  c(1, sqrt2plus1, -1, -1),
  c(-1, -sqrt2plus1, 1, -1),
  c(1, -sqrt2plus1, 1, -1),
  c(-1, sqrt2plus1, 1, -1),
  c(1, sqrt2plus1, 1, -1),
  c(-1, -sqrt2plus1, -1, 1),
  c(1, -sqrt2plus1, -1, 1),
  c(-1, sqrt2plus1, -1, 1),
  c(1, sqrt2plus1, -1, 1),
  c(-1, -sqrt2plus1, 1, 1),
  c(1, -sqrt2plus1, 1, 1),
  c(-1, sqrt2plus1, 1, 1),
  c(1, sqrt2plus1, 1, 1),
  c(-sqrt2plus1, -1, -1, -1),
  c(sqrt2plus1, -1, -1, -1),
  c(-sqrt2plus1, 1, -1, -1),
  c(sqrt2plus1, 1, -1, -1),
  c(-sqrt2plus1, -1, 1, -1),
  c(sqrt2plus1, -1, 1, -1),
  c(-sqrt2plus1, 1, 1, -1),
  c(sqrt2plus1, 1, 1, -1),
  c(-sqrt2plus1, -1, -1, 1),
  c(sqrt2plus1, -1, -1, 1),
  c(-sqrt2plus1, 1, -1, 1),
  c(sqrt2plus1, 1, -1, 1),
  c(-sqrt2plus1, -1, 1, 1),
  c(sqrt2plus1, -1, 1, 1),
  c(-sqrt2plus1, 1, 1, 1),
  c(sqrt2plus1, 1, 1, 1)
)

tetrahedra =  1 + c(
  c(38, 22, 6, 54),
  c(30, 14, 46, 62),
  c(63, 31, 47, 15),
  c(39, 55, 7, 23),
  c(26, 10, 42, 58),
  c(32, 16, 0, 48),
  c(19, 35, 51, 3),
  c(44, 12, 28, 60),
  c(27, 59, 43, 11),
  c(18, 2, 50, 34),
  c(45, 29, 61, 13),
  c(33, 17, 1, 49),
  c(41, 25, 9, 57),
  c(37, 5, 53, 21),
  c(40, 24, 56, 8),
  c(36, 4, 52, 20)
)

pts4 <- vertices[tetrahedra, ]

stereog <- function(q, r2 = 3 + sqrt2plus1^2){
  acos(q[4]/sqrt(r2))/pi * q[1:3]/sqrt(r2-q[4]*q[4])
}

pts <- t(apply(pts4, 1L, stereog))

d <- delaunay(pts)
cxfacets <- cxhull::cxhull(pts)$facets

edges <- d$exteriorEdges

edges3 <- list()

j <- 0
for(edge in edges){
  j <- j + 1
  A <- edge$A
  x <- sapply(cxfacets, function(f){
    c(crossprod(f$normal, A)) + f$offset
  })
  Abelongs <- which(abs(x) < sqrt(.Machine$double.eps))
  B <- edge$B
  if(length(Abelongs)){
    x <- sapply(Abelongs, function(i){
      f <- cxfacets[[i]]
      c(crossprod(f$normal, B)) + f$offset
    })
    Bext <- any(abs(x) < sqrt(.Machine$double.eps))
    if(Bext){
      edges3 <- append(edges3, list(edge))
    }
  }
}




library(gMOIP)
hull <- geometry::convhulln(pts)
edges3 <- list()
for(edge in edges){
  if(inHull(rbind(edge$A), pts, hull, tol = sqrt(.Machine$double.eps)) == 0 &&
     inHull(rbind(edge$B), pts, hull, tol = sqrt(.Machine$double.eps)) == 0){
    edges3 <- append(edges3, list(edge))
  }
}
plotDelaunay3D(d, exteriorEdgesAsTubes = TRUE, tubeRadius = 0.03, tubeColor = "green")
invisible(lapply(edges3, function(edge) edge$plot(edgeAsTube = TRUE, tubeRadius = 0.03, tubeColor = "green")))




library(geometry)
ch <- convhulln(pts)

M <- lapply(edges, function(edge){
  rbind(edge$A, edge$B)
})


edges3 <- list()
for(edge in edges){
  if(inhulln(ch, rbind(edge$A)) && inhulln(ch, rbind(edge$B))){
    edges3 <- append(edges3, list(edge))
  }
}


plotDelaunay3D(d, exteriorEdgesAsTubes = TRUE, tubeRadius = 0.03, tubeColor = "orange")









library(cxhull)
h <- cxhull(pts, triangulate = TRUE)
library(rgl)
for(i in 1:nrow(h$edges)){
  edge <- h$edges[i,]
  tube <- cylinder3d(rbind(pts[edge[1],], pts[edge[2],]), sides= 30, radius = 0.035)
  shade3d(tube, color= "green")
}
