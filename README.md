# spatstat.geom

## Spatial data classes and geometrical functionality for the spatstat family

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/spatstat.geom)](http://CRAN.R-project.org/package=spatstat.geom) 
[![GitHub R package version](https://img.shields.io/github/r-package/v/spatstat/spatstat.geom)](https://github.com/spatstat/spatstat.geom)

The original `spatstat` package has been split into
several sub-packages. This package `spatstat.geom` is one of the sub-packages.
It defines classes of spatial objects, and supports
geometrical operations on them.

### Overview

`spatstat.geom` supports

- classes of spatial objects
- basic support for spatial objects
- graphics
- interaction with spatial objects
- geometrical operations
- morphological operations
- image processing

### Detailed contents

For a full list of functions, see the help file for `spatstat.geom-package`.

#### Classes of spatial objects in 2D

- point patterns (`ppp`)
- spatial regions or windows (`owin`)
- pixel images (`im`)
- line segment patterns (`psp`)
- tessellations (`tess`)
- linear networks (`linnet`)

#### Classes of spatial objects in 3D and higher

- 3D point patterns (`pp3`)
- general multidimensional space-time point patterns (`ppx`)

#### Basic support

- printing spatial objects
- basic summary information about spatial objects
- subsetting, splitting, superimposing spatial patterns
- extracting information such as spatial coordinates
- pixellation
- discretisation of coordinates
- interpolation of images
- solution set, level set
- identifying duplicated elements

#### Graphics

- plotting spatial patterns
- plotting images (depicted as colour images, contour plots, perspective views, composite images)
- graphical information (layered objects, colour maps, symbol maps)

#### Interactive

- interactive data entry and editing for spatial objects
- `identify` methods for identifying individual elements
- simple interactive panels (`simplepanel`)

#### Geometry

- geometrical transformations (`rotate`, `scalardilation`,
`shift`, `reflect`, `flipxy`, `affine`)
- set operations (`intersect.owin`, `union.owin`, `complement.owin`, `setminus.owin`)
- test whether a point falls inside a set
- convex hull
- geometrical mensuration (measuring areas, lengths, angles, diameter)

#### Distance Operations

- nearest neighbour
- nearest neighbour distance
- pairwise distances
- distance map/ distance transform

#### Geometrical constructions

- tessellations (`tess`, `hextess`, `quadrats`, `venn.tess`, `polartess`, `dirichlet`, `delaunay`, `quantess`, `rpoislinetess`)
- triangulation
- discretisation
- bounding box

#### Mathematical Morphology

- erosion, dilation
- opening, closing
- morphological distance transform

#### Image Processing

- pixel-by-pixel arithmetic
- set covariance function, convolution of images

#### Elementary summary statistics

- quadrat counting

The reorganisation of `spatstat` into a family of packages is described
on the GitHub repository
[spatstat/spatstat](https://github.com/spatstat/spatstat).

This repository holds the latest development version of `spatstat.geom`.
For the latest public release on CRAN, click the green badge above.
