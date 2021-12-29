library(tessellation)

vectorcycle <- function(v, i){
  c(tail(v, i), head(v, length(v)-i))
}

alpha_ <- seq(0, 2*pi, length.out = 361)[-1L]
circle <- cbind(14.8*cos(alpha_), 14.8*sin(alpha_))

library(paletteer)
colors <- paletteer_c("pals::kovesi.cyclic_mygbm_30_95_c78_s25", 360)#281)

theta <- seq(0, 110, length.out = 382L) # 100 300
x <- sqrt(theta) * cos(theta)
y <- sqrt(theta) * sin(theta)
pts <- cbind(x,y)
opar <- par(mar = c(0, 0, 0, 0), bg = "black")
# Here is a Fermat spiral:
#plot(pts, asp = 1, xlab = NA, ylab = NA, axes = FALSE, pch = 19, col = "white")
# And here is its Voronoï diagram:
svg("zz.svg")
opar <- par(mar = c(0, 0, 0, 0), bg = "transparent")
plot(NULL, asp = 1, xlim = c(-15, 15), ylim = c(-15, 15),
     xlab = NA, ylab = NA, axes = FALSE)
del <- delaunay(pts)
v <- voronoi(del)
attr(v, "nbounded") # 281
polygon(circle, col = "black")
plotVoronoiDiagram(v, colors = vectorcycle(colors, 100L))
dev.off()


Fplot <- function(i){
  opar <- par(mar = c(0, 0, 0, 0), bg = "transparent")
  plot(NULL, asp = 1, xlim = c(-15, 15), ylim = c(-15, 15),
       xlab = NA, ylab = NA, axes = FALSE)
  polygon(circle, col = "black")
  plotVoronoiDiagram(v , colors = vectorcycle(colors, i))
}

library(rsvg)
for(i in 1L:360L){
  svgname <- sprintf("vfermat%03d.svg", i)
  svg(svgname, bg = "transparent")
  Fplot(i)
  dev.off()
  pngname <- sprintf("vfermat%03d.png", i)
  rsvg_png(svgname, file = pngname, width = 512, height = 512)
  system(sprintf("mogrify -rotate %d %s", i, pngname))
}

for(i in 1L:360L){
  svgname <- sprintf("vfermat%03d.svg", i)
  pngname <- sprintf("vvfermat%03d.png", i)
  rsvg_png(svgname, file = pngname, width = 512, height = 512)
  system(sprintf("mogrify -distort SRT %d %s", i, pngname))
}



pngfiles <- list.files(pattern = "^vvfermat.*png$")

library(gifski)
gifski(
  pngfiles,
  gif_file = "VoronoiFermatSpiral3.gif",
  width = 512,
  height = 512,
  delay = 1/30
)

