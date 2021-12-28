# tessellation

<!-- badges: start -->
[![R-CMD-check](https://github.com/stla/tessellation/workflows/R-CMD-check/badge.svg)](https://github.com/stla/tessellation/actions)
<!-- badges: end -->

This package is not on CRAN yet. You can install it by typing:

```r
remotes::install_github("stla/tessellation", dependencies = TRUE, build_vignettes = TRUE)
```

## Delaunay and Voronoï tessellations.

The Voronoï cell of a point inside the Utah teapot:

![](https://raw.githubusercontent.com/stla/tessellation/main/inst/screenshots/UtahTeapot.png)

The Voronoï diagram (restricted to its bounded cells) of a cube surrounded by three perpendicular circles:

![](https://raw.githubusercontent.com/stla/tessellation/main/inst/screenshots/surroundedCube.png)

The Voronoï diagram (restricted to its bounded cells) of a dodecahedron surrounded by three perpendicular circles:

![](https://raw.githubusercontent.com/stla/tessellation/main/inst/screenshots/dodecahedron.gif)

A strange Voronoï cell:

![](https://raw.githubusercontent.com/stla/tessellation/main/inst/screenshots/strangeVoronoiCell.gif)

In fact this is not a Voronoï cell, I did a mistake :). Voronoï cells are convex.
