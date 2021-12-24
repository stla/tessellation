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

plotDelaunay3D(d, exteriorEdgesAsTubes = TRUE, tubeRadius = 0.03, tubeColor = "orange")
