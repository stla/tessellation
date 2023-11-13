# tessellation 2.3.0

New function `cellVolume`. It returns the volume of a bounded 2D or 3D Voronoï
cell (i.e. the area in 2D, the volume in 3D). 
For a bounded 3D Voronoï cell, in addition to the volume of the cell, the 
surface area of the cell is also returned.


# tessellation 2.2.0

- Fixed a misspelled word.

- The package does no longer depend on the 'interp' package.

- The package does no longer depend on the 'randomcoloR' package, but it now 
depends on the packages 'colorsGen' and 'Polychrome'.


# tessellation 2.1.3

Fixed a minor mistake in the `plotVoronoiDiagram` function.


# tessellation 2.1.2

Fixed a small mistake in the C code.


# tessellation 2.1.1

Changed `sprintf` to `snprintf` in the C code.


# tessellation 2.1.0

* Improved doc and examples.

* Now the elevated Delaunay tessellation also returns the total area of the 
triangulation.


# tessellation 2.0.0

New feature: elevated Delaunay tessellation (also called 2.5D Delaunay 
triangulation).


# tessellation 1.0.0

First release.
