

theta = seq(0, 100, length.out = 300)
x = sqrt(theta) * cos(theta)
y = sqrt(theta) * sin(theta)

pts = cbind(x,y)

opar <- par(mar = c(0, 0, 0, 0), bg = "black")
plot(NULL, asp = 1, xlim = c(-15, 15), ylim = c(-15, 15),
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
