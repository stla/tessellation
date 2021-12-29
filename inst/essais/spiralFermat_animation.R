library(tessellation)

vectorcycle <- function(v, i){
  c(tail(v, i), head(v, length(v)-i))
}

library(paletteer)
colors <- paletteer_c("pals::kovesi.cyclic_mygbm_30_95_c78_s25", 360)#281)

theta <- seq(0, 110, length.out = 382L) # 100 300
x <- sqrt(theta) * cos(theta)
y <- sqrt(theta) * sin(theta)
pts <- cbind(x,y)
opar <- par(mar = c(0, 0, 0, 0), bg = "black")
# Here is a Fermat spiral:
plot(pts, asp = 1, xlab = NA, ylab = NA, axes = FALSE, pch = 19, col = "white")
# And here is its Voronoï diagram:
svg("zz.svg")
opar <- par(mar = c(0, 0, 0, 0), bg = "black")
plot(NULL, asp = 1, xlim = c(-15, 15), ylim = c(-15, 15),
     xlab = NA, ylab = NA, axes = FALSE)
del <- delaunay(pts)
v <- voronoi(del)
attr(v, "nbounded") # 281
plotVoronoiDiagram(v, colors = vectorcycle(colors, 100L))
dev.off()


Fplot <- function(i){
  opar <- par(mar = c(0, 0, 0, 0), bg = "black")
  plot(NULL, asp = 1, xlim = c(-15, 15), ylim = c(-15, 15),
       xlab = NA, ylab = NA, axes = FALSE)
  plotVoronoiDiagram(v , colors = vectorcycle(colors, i))
}

library(tikzDevice)
library(stringr)
for(i in 1L:281L){
  fplot <- function() Fplot(i)
  plot2tikz(
    fplot, compile=FALSE,
    filename = sprintf("fermat%03d", i),
    documentDeclaration = "\\documentclass[12pt]{standalone}\n",
    width=6, height=6,
    bg="black", fg="black"
  )
}

texfiles <- list.files(pattern = "^fermat")
for(texfile in texfiles){
  tools::texi2dvi(texfile)
}

dvifiles <- list.files(pattern = "^fermat.*dvi$")
for(dvifile in dvifiles){
  command <- sprintf("dvipng -T 512px,512px %s", dvifile)
  system(command)
}

system("mogrify -resize 512x512! fermat*.png")


pngfiles <- list.files(pattern = "^fermat.*png$")

library(gifski)
gifski(
  pngfiles,
  gif_file = "VoronoiFermatSpiral.gif",
  width = 512,
  height = 512,
  delay = 1/30
)

