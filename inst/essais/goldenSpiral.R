
phi <- (1+sqrt(5)) / 2

theta = seq(1, 30, length.out = 100)
x = cos(theta) * phi^(theta/pi)
y = sin(theta) * phi^(theta/pi)

pts = cbind(x,y)

opar <- par(mar = c(0, 0, 0, 0), bg = "black")
plot(pts, asp = 1, # xlim = c(-15, 15), ylim = c(-15, 15),
     xlab = NA, ylab = NA, axes = FALSE)

del <- delaunay(pts)
v <- voronoi(del)
plotVoronoiDiagram(v , colors = viridisLite::turbo(281))
plotVoronoiDiagram(v , colors = randomcoloR::distinctColorPalette(281))
plotVoronoiDiagram(v, luminosity = "dark")


fplot <- function(){
  opar <- par(mar = c(0, 0, 0, 0), bg = "black")
  plot(NULL, asp = 1, xlim = c(-15, 15), ylim = c(-15, 15),
       xlab = NA, ylab = NA, axes = FALSE)
  plotVoronoiDiagram(v , colors = viridisLite::turbo(281))
}
