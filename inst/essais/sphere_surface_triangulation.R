library(tessellation)

n <- 10
rho <- 1
theta <- seq(0, 2*pi, length.out = n) # azimuthal coordinate running from 0 to 2*pi
phi <- seq(0, pi, length.out = n) # polar coordinate running from 0 to pi (colatitude)
grd <- expand.grid(theta=theta, phi=phi)

x <- rho * cos(grd$theta) * sin(grd$phi)
y <- rho * sin(grd$theta) * sin(grd$phi)
z <- rho * cos(grd$phi)
xyz <- cbind(x,y,z)


library(cxhull)
hull <- cxhull(xyz, triangulate = TRUE)

library(rgl)
open3d()
for(facet in hull$facets){
  triangles3d(xyz[facet[["vertices"]], ], color = "green")
}
