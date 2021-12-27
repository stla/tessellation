library(CVXR)
library(rgl)
library(cxhull)
cube <- t(cube3d()$vb[-4,])
pts <- rbind(cube, c(0,1,1))

