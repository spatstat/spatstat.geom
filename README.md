# spatstat.geom

## Spatial data classes and geometrical functionality for the spatstat family

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/spatstat.geom)](http://CRAN.R-project.org/package=spatstat.geom) 
[![GitHub R package version](https://img.shields.io/github/r-package/v/spatstat/spatstat.geom)](https://github.com/spatstat/spatstat.geom)


You are viewing the GitHub repository which
holds the latest **development version** of `spatstat.geom`.
For the latest public release on CRAN, click the green badge above.

Contents:

 - [Overview of `spatstat.geom`](#overview)
 - [Detailed contents of package](#detailed)
 - [Installing the package](#installing)
 - [Bug reports](#bugreports)
 - [Questions](#questions)
 - [Proposing changes to code](#proposing)
 - [Future development](#future)

___

### <a name="overview"></a> Overview of `spatstat.geom`

The original `spatstat` package has been split into
several sub-packages, listed on the GitHub repository
[spatstat/spatstat](https://github.com/spatstat/spatstat).

This package `spatstat.geom` is one of the sub-packages.
It defines **classes of spatial objects**, and supports
**geometrical operations** on them. It provides

- classes of spatial objects
- basic support for spatial objects
- graphics
- interaction with spatial objects
- geometrical operations
- morphological operations
- image processing

___

### <a name="detailed"></a> Detailed contents of package

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

- nearest neighbour distance (`nndist`)
- find the nearest neighbour (`nnwhich`)
- pairwise distances (`pairdist`)
- nearest neighbour from one pattern to another (`nncross`)
- pairwise distances between one pattern and another (`crossdist`)
- distance transform (`distmap`, `distfun`)
- non-Euclidean distance (`convexmetric`)

#### Geometrical constructions

- creating tessellations (`tess`, `hextess`, `quadrats`, `venn.tess`, `polartess`, `dirichlet`, `delaunay`, `quantess`)
- triangulation
- discretisation
- bounding box

#### Mathematical Morphology

- erosion, dilation
- opening, closing
- morphological distance transform

#### Image Processing

- pixel-by-pixel computation (`Math.im`, `eval.im`, `im.apply`)
- set covariance function (`setcov`)
- convolution of images (`imcov`)

#### Elementary summary statistics

Summary statistics are mostly provided in the `spatstat.explore` package.
However, `spatstat.geom` provides functions for calculating

- average intensity of a point pattern (`intensity`)
- quadrat counts (`quadratcount`)

#### Elementary generation of random patterns

Functions for generating random spatial patterns
are mostly provided in the `spatstat.random` package.
However, `spatstat.geom` provides basic functions for generating

- random points in a rectangle (`runifrect`)
- quasirandom points (`rQuasi`)

___

### <a name="installing"></a> Installing the package

This repository contains the _development version_ of
`spatstat.geom`. The easiest way to install the development version
is to start R and type

```R
repo <- c('https://spatstat.r-universe.dev', 'https://cloud.r-project.org')
install.packages("spatstat.geom", dependencies=TRUE, repos=repo)
```

To install the latest _public release_ of `spatstat.geom`,
type

```R
install.packages("spatstat.geom")
```

___

## <a name="bugreports"></a> Bug reports 

Users are encouraged to report bugs.
If you find a bug in a `spatstat` function,
please identify the sub-package containing that function.
Visit the GitHub repository for the sub-package, 
click the `Issues` tab at the top of the page, 
and press *new issue* to start a new bug report, documentation correction
or feature request.

**Please do not post questions** on the Issues pages,
because they are too clunky for correspondence.

## <a name="questions"></a> Questions about spatstat

For questions about the `spatstat` package family, first check 
the question-and-answer website
[stackoverflow](http://stackoverflow.com/questions/tagged/spatstat)
to see whether your question has already been asked and answered.
If not, you can either post your question at stackoverflow, or
email the authors.

## <a name="proposing"></a> Proposing changes to the code

Feel free to fork `spatstat.geom`, make changes to the code,
and ask us to include them in the package by making a github *pull request*. 

## <a name="future"></a> Future development

The `spatstat` package family is the result of 30 years of software development
and contains over 200,000 lines of code.
It is still under development,
motivated by the needs of researchers in many fields,
and driven by innovations in statistical science.
We welcome contributions of code, and suggestions
for improvements.
