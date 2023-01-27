# tessellation

<!-- badges: start -->
[![R-CMD-check](https://github.com/stla/tessellation/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/stla/tessellation/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

To install the development version, run:

```r
remotes::install_github("stla/tessellation", dependencies = TRUE)
```

I did a vignette but it was too big for the CRAN submission. So I used 
**pkgdown**, and you can find the [vignette here](https://stla.github.io/tessellation/articles/the-tessellation-package.html).


## Delaunay and Voronoï tessellations.

The Voronoï cell of a point inside the Utah teapot:

![](https://raw.githubusercontent.com/stla/tessellation/main/inst/screenshots/UtahTeapot.png)

The Voronoï diagram (restricted to its bounded cells) of a cube surrounded by three perpendicular circles:

![](https://raw.githubusercontent.com/stla/tessellation/main/inst/screenshots/surroundedCube.png)

The Voronoï diagram (restricted to its bounded cells) of a dodecahedron surrounded by three perpendicular circles:

![](https://raw.githubusercontent.com/stla/tessellation/main/inst/screenshots/dodecahedron.gif)

A strange Voronoï cell:

![](https://raw.githubusercontent.com/stla/tessellation/main/inst/screenshots/strangeVoronoiCell.gif)

The Voronoï diagram of some points along a Fermat spiral:

![](https://raw.githubusercontent.com/stla/tessellation/main/inst/screenshots/VoronoiFermatSpiral.gif)

![](https://raw.githubusercontent.com/stla/tessellation/main/inst/screenshots/VoronoiFermatSpiral2.gif)

![](https://raw.githubusercontent.com/stla/tessellation/main/inst/screenshots/VoronoiFermatSpiral3.gif)

An elevated Delaunay tessellation:

![](https://raw.githubusercontent.com/stla/tessellation/main/inst/screenshots/volcano.png)
