\name{spatstat.geom-package}
\alias{spatstat.geom-package} 
\alias{spatstat.geom} 
\docType{package}
\title{The spatstat.geom Package}
\description{
  The \pkg{spatstat.geom} package
  belongs to the \pkg{spatstat} family of packages.
  It defines classes of geometrical objects such as windows
  and point patterns, and provides functionality for geometrical
  operations on them.
}
\details{
  \pkg{spatstat} is
  a family of \R packages
  for the statistical analysis of spatial data.
  Its main focus is the analysis of
  spatial patterns of points in two-dimensional space.

  The original \pkg{spatstat} package
  has now been split into several
  sub-packages.

  
  This sub-package \pkg{spatstat.geom} defines the main classes
  of geometrical objects (such as windows, point patterns, 
  line segment patterns, pixel images) and supports geometrical
  operations (such as shifting and rotating, measuring areas and
  distances, finding nearest neighbours in a point pattern).

  Functions for performing statistical analysis and modelling
  are in the separate sub-packages \pkg{spatstat.explore} and \pkg{spatstat.model}.
  
  Functions for linear networks
  are in the separate sub-package \pkg{spatstat.linnet}.

  For an overview of all the functions available in the \pkg{spatstat} family,
  see the help file for \pkg{spatstat} in the \pkg{spatstat} package.
}
\section{Structure of the spatstat family}{  
  The original \pkg{spatstat} package grew to be very large,
  and CRAN requested that the package be divided 
  into several \bold{sub-packages}.
  Currently the sub-packages are:
  \itemize{
    \item \pkg{spatstat.utils} containing basic utilities
    \item \pkg{spatstat.data} containing datasets
    \item \pkg{spatstat.sparse} containing linear algebra utilities
    \item \pkg{spatstat.geom} containing geometrical objects
    and geometrical operations
    \item \pkg{spatstat.random} containing code for generating
    random spatial patterns
    \item \pkg{spatstat.explore} containing the main functionality
    for exploratory and non-parametric analysis of spatial data
    \item \pkg{spatstat.model} containing the main functionality
    for statistical modelling and inference for spatial data
    \item \pkg{spatstat.linnet} containing functions for
    spatial data on a linear network
    \item \pkg{spatstat}, which simply loads the other sub-packages
    listed above, and provides documentation.
  }
  When you install \pkg{spatstat}, these sub-packages are also
  installed. Then if you load the \pkg{spatstat} package by typing
  \code{library(spatstat)}, the other sub-packages listed above will
  automatically be loaded or imported.
  For an overview of all the functions available in these sub-packages,
  see the help file for \pkg{spatstat} in the \pkg{spatstat} package,

  Additionally there are several \bold{extension packages:}
  \itemize{
    \item \pkg{spatstat.gui} for interactive graphics
    \item \pkg{spatstat.local} for local likelihood
    (including geographically weighted regression)
    \item \pkg{spatstat.Knet} for additional, computationally efficient code
    for linear networks
    \item \pkg{spatstat.sphere} (under development) for spatial data
    on a sphere, including spatial data on the earth's surface
  }
  The extension packages must be installed separately
  and loaded explicitly if needed. They also have separate documentation.
}
\section{OVERVIEW OF CAPABILITIES}{
  Following is an overview of the capabilities of
  the \pkg{spatstat.geom} sub-package.

  \bold{Types of spatial data:}

  The main types of spatial data supported by \pkg{spatstat.geom} are:

  \tabular{ll}{
    \code{\link{ppp}} \tab point pattern \cr
    \code{\link{owin}} \tab window (spatial region) \cr
    \code{\link{im}} \tab pixel image \cr
    \code{\link{psp}} \tab line segment pattern \cr
    \code{\link{tess}} \tab tessellation \cr
    \code{\link{pp3}} \tab three-dimensional point pattern \cr
    \code{\link{ppx}} \tab point pattern in any number of dimensions \cr
  }
  Additional data types are supported in \pkg{spatstat.linnet}.
  
  \bold{To create a point pattern:}
  
  \tabular{ll}{
    \code{\link{ppp}} \tab
    create a point pattern from \eqn{(x,y)} and window information
    \cr
    \tab
    \code{ppp(x, y, xlim, ylim)} for rectangular window\cr
    \tab
    \code{ppp(x, y, poly)} for polygonal window \cr
    \tab
    \code{ppp(x, y, mask)} for binary image window \cr
    \code{\link{as.ppp}} \tab
    convert other types of data to a \code{ppp} object \cr
    \code{\link{clickppp}} \tab
    interactively add points to a plot \cr
    \code{\link{marks<-}}, \code{\%mark\%}  \tab
    attach/reassign marks to a point pattern
  }
      
  \bold{To simulate a random point pattern:}

  Most of the methods for generating random data are
  provided in \pkg{spatstat.random}. The following basic methods
  are supplied in \pkg{spatstat.geom}:
  
  \tabular{ll}{
	\code{\link{runifrect}} \tab
	generate \eqn{n} independent uniform random points in a	rectangle \cr
	\code{\link{rsyst}} \tab
	systematic random sample of points \cr
        \code{\link{rjitter}} \tab
	apply random displacements to points in a pattern\cr
      }

   \bold{Standard point pattern datasets:}

   Datasets installed in the \pkg{spatstat} family are provided
   in the sub-package \code{spatstat.data}. 

   \bold{To manipulate a point pattern:}

      \tabular{ll}{
	\code{\link{plot.ppp}} \tab
	plot a point pattern (e.g. \code{plot(X)}) \cr
	\code{spatstat.gui::iplot} \tab
	plot a point pattern interactively \cr 
	\code{\link{edit.ppp}} \tab interactive text editor \cr
	\code{\link{[.ppp}} \tab 
	extract or replace a subset of a point pattern \cr
        \tab \code{pp[subset]} or \code{pp[subwindow]} \cr
	\code{\link{subset.ppp}} \tab
	extract subset of point pattern satisfying a condition \cr
	\code{\link{superimpose}} \tab
	combine several point patterns  \cr
	\code{\link{by.ppp}} \tab
	apply a function to sub-patterns of a point pattern \cr
	\code{\link{cut.ppp}} \tab
	classify the points in a point pattern \cr
	\code{\link{split.ppp}} \tab
	divide pattern into sub-patterns \cr 
	\code{\link{unmark}} \tab
	remove marks  \cr
	\code{\link{npoints}} \tab
	count the number of points  \cr
	\code{\link{coords}} \tab
	extract coordinates, change coordinates  \cr
	\code{\link{marks}} \tab
	extract marks, change marks or attach marks  \cr
	\code{\link{rotate}} \tab
	rotate pattern  \cr
	\code{\link{shift} } \tab
	translate pattern  \cr
	\code{\link{flipxy} } \tab
	swap \eqn{x} and \eqn{y} coordinates  \cr
	\code{\link{reflect} } \tab
	reflect in the origin  \cr
	\code{\link{periodify} } \tab
	make several translated copies  \cr
	\code{\link{affine}} \tab
	apply affine transformation\cr
	\code{\link{scalardilate}} \tab
	apply scalar dilation\cr
	\code{\link{nnmark}} \tab
	mark value of nearest data point\cr
	\code{\link{identify.ppp}} \tab
	interactively identify points \cr
	\code{\link{unique.ppp}} \tab
	remove duplicate points \cr
	\code{\link{duplicated.ppp}} \tab
	determine which points are duplicates \cr
	\code{\link{uniquemap.ppp}} \tab
	map duplicated points to unique points \cr
	\code{\link{connected.ppp}} \tab find clumps of points \cr
	\code{\link{dirichlet}} \tab
	compute Dirichlet-Voronoi tessellation \cr
	\code{\link{delaunay}} \tab
	compute Delaunay triangulation \cr 
	\code{\link{delaunayDistance}} \tab
	graph distance in Delaunay triangulation \cr 
	\code{\link{convexhull}} \tab compute convex hull \cr
	\code{\link{discretise}} \tab discretise coordinates \cr
	\code{\link{pixellate.ppp}} \tab approximate point pattern by 
	pixel image \cr
	\code{\link{as.im.ppp}} \tab approximate point pattern by 
	pixel image 
      }
      See \code{\link{spatstat.options}} to control plotting behaviour.
      
    \bold{To create a window:}

    An object of class \code{"owin"} describes a spatial region
      (a window of observation).

      \tabular{ll}{
	\code{\link{owin}}	\tab	Create a window object \cr
	\tab \code{owin(xlim, ylim)} for rectangular window \cr
	\tab \code{owin(poly)} for polygonal window \cr
	\tab \code{owin(mask)} for binary image window \cr
	\code{\link{Window}}	\tab
	Extract window of another object \cr
	\code{\link{Frame}}	\tab
	Extract the containing rectangle ('frame') of another object \cr
	\code{\link{as.owin}}	\tab
	Convert other data to a window object \cr
	\code{\link{square}}    \tab	make a square window \cr
	\code{\link{disc}}    \tab	make a circular window \cr
	\code{\link{ellipse}}    \tab	make an elliptical window \cr
	\code{\link{ripras}}    \tab
	Ripley-Rasson estimator of window, given only the points \cr
	\code{\link{convexhull}} \tab compute convex hull of something \cr 
	\code{\link[spatstat.data]{letterR}}    \tab
	polygonal window in the shape of the \R logo \cr
	\code{\link{clickpoly}}    \tab
	interactively draw a polygonal window  \cr
	\code{\link{clickbox}}    \tab
	interactively draw a rectangle  
      }

    \bold{To manipulate a window:}

    \tabular{ll}{
	\code{\link{plot.owin}}	\tab	plot a window. \cr
	\tab		\code{plot(W)}\cr
	\code{\link{boundingbox}} \tab
	Find a tight bounding box for the window \cr
	\code{\link{erosion}}	\tab
	erode window by a distance r\cr
	\code{\link{dilation}}	\tab
	dilate window by a distance r\cr
	\code{\link{closing}}	\tab
	close window by a distance r\cr
	\code{\link{opening}}	\tab
	open window by a distance r\cr
	\code{\link{border}}	\tab
	difference between window and its erosion/dilation \cr 
	\code{\link{complement.owin}}	\tab
	invert (swap inside and outside)\cr
	\code{\link{simplify.owin}}	\tab
	approximate a window by a simple polygon  \cr
	\code{\link{rotate}} \tab rotate window  \cr
	\code{\link{flipxy}} \tab swap \eqn{x} and \eqn{y} coordinates  \cr
	\code{\link{shift} } \tab translate window  \cr
	\code{\link{periodify} } \tab make several translated copies  \cr
	\code{\link{affine}} \tab apply affine transformation \cr
	\code{\link{as.data.frame.owin}} \tab
	convert window to data frame 
      }

    \bold{Digital approximations:}

    \tabular{ll}{
	\code{\link[spatstat.geom]{as.mask}}	\tab
	Make a discrete pixel approximation of a given window \cr
	\code{\link{as.im.owin}} \tab convert window to pixel image \cr
	\code{\link{pixellate.owin}} \tab convert window to pixel image \cr
	\code{\link{commonGrid}} \tab find common pixel grid for windows \cr
	\code{\link{nearest.raster.point}} \tab
	map continuous coordinates to raster locations\cr
	\code{\link{raster.x}} \tab
	raster x coordinates \cr
	\code{\link{raster.y}} \tab
	raster y coordinates \cr
	\code{\link{raster.xy}} \tab
	raster x and y coordinates \cr
	\code{\link{as.polygonal}} \tab
	convert pixel mask to polygonal window 
      }
      See \code{\link{spatstat.options}} to control the approximation

    \bold{Geometrical computations with windows:}

    \tabular{ll}{
	\code{\link{edges}}	\tab	extract boundary edges \cr
	\code{\link{intersect.owin}}	\tab	intersection of two windows\cr
	\code{\link{union.owin}}	\tab	union of two windows\cr
	\code{\link{setminus.owin}}	\tab	set subtraction of two windows\cr
	\code{\link{inside.owin}}	\tab	determine whether a point is inside a window\cr
	\code{\link{area.owin}}	\tab	compute area \cr
	\code{\link{perimeter}}	\tab	compute perimeter length \cr
	\code{\link{diameter.owin}}	\tab	compute diameter\cr
	\code{\link{incircle}}	\tab	find largest circle inside a window \cr
	\code{\link{inradius}}	\tab	radius of incircle \cr
	\code{\link{connected.owin}}    \tab find connected components of window \cr
	\code{\link{eroded.areas}}	\tab	compute areas of eroded windows\cr
	\code{\link{dilated.areas}}	\tab	compute areas of dilated windows\cr
	\code{\link{bdist.points}}	\tab	compute distances from data points to window boundary \cr
	\code{\link{bdist.pixels}}	\tab	compute distances from all pixels to window boundary \cr
	\code{\link{bdist.tiles}}	\tab
	boundary distance for each tile in tessellation \cr
	\code{\link{distmap.owin}}	\tab	distance transform image \cr
	\code{\link{distfun.owin}}	\tab	distance transform \cr
	\code{\link{centroid.owin}}	\tab compute centroid (centre of mass) of window\cr
	\code{\link{is.subset.owin}}    \tab determine whether one
	window contains another \cr
	\code{\link{is.convex}} \tab determine whether a window is convex \cr 
	\code{\link{convexhull}} \tab compute convex hull \cr
	\code{\link{triangulate.owin}} \tab decompose into triangles \cr
	\code{\link[spatstat.geom]{as.mask}} \tab pixel approximation of window \cr
	\code{\link{as.polygonal}} \tab polygonal approximation of window \cr
	\code{\link{is.rectangle}} \tab test whether window is a rectangle \cr 
	\code{\link{is.polygonal}} \tab test whether window is polygonal \cr 
	\code{\link{is.mask}} \tab test whether window is a mask \cr 
	\code{\link{setcov}} \tab spatial covariance function of window \cr
	\code{\link{pixelcentres}} \tab extract centres of pixels in mask \cr
	\code{\link{clickdist}}    \tab
	measure distance between two points clicked by user
      }

    \bold{Pixel images:}
    An object of class \code{"im"} represents a pixel image. 
    Such objects are returned by some of the functions in
      \pkg{spatstat} including \code{\link[spatstat.explore]{Kmeasure}},
      \code{\link{setcov}} and \code{\link[spatstat.explore]{density.ppp}}. 
      \tabular{ll}{
	\code{\link{im}} \tab create a pixel image\cr
	\code{\link{as.im}} \tab convert other data to a pixel image\cr
	\code{\link{pixellate}} \tab convert other data to a pixel image\cr
	\code{\link{as.matrix.im}} \tab convert pixel image to matrix\cr
	\code{\link{as.data.frame.im}} \tab convert pixel image to data frame\cr
	\code{\link{as.function.im}} \tab convert pixel image to function\cr
	\code{\link{plot.im}}	\tab	plot a pixel image on screen as a digital image\cr
	\code{\link{contour.im}}	\tab draw contours of a pixel image \cr
	\code{\link{persp.im}}	\tab draw perspective plot of a pixel image \cr
	\code{\link{rgbim}}	\tab create colour-valued pixel image \cr
	\code{\link{hsvim}}	\tab create colour-valued pixel image \cr
	\code{\link{[.im}} 	\tab extract a subset of a pixel image\cr
	\code{\link{[<-.im}} 	\tab replace a subset of a pixel image\cr
	\code{\link{rotate.im}} \tab rotate pixel image \cr
	\code{\link{shift.im}} \tab apply vector shift to pixel image \cr
	\code{\link{affine.im}} \tab apply affine transformation to image \cr
	\code{X}	\tab	print very basic information about image \code{X}\cr
	\code{\link{summary}(X)} \tab	summary of image \code{X} \cr
	\code{\link{hist.im}} \tab	histogram of image \cr
	\code{\link{mean.im}} \tab	mean pixel value of image  \cr
	\code{\link{integral.im}} \tab	integral of pixel values  \cr
	\code{\link{quantile.im}} \tab	quantiles of image  \cr
	\code{\link{cut.im}} \tab	convert numeric image to factor image \cr
	\code{\link{is.im}} \tab test whether an object is a pixel image\cr
	\code{\link{interp.im}} \tab interpolate a pixel image\cr
	\code{\link{connected.im}} \tab find connected components \cr
	\code{\link{compatible.im}} \tab test whether two images have
	compatible dimensions \cr
	\code{\link{harmonise.im}} \tab make images compatible \cr
	\code{\link{commonGrid}} \tab find a common pixel grid for images \cr
	\code{\link{eval.im}} \tab evaluate any expression involving images\cr
	\code{\link{im.apply}} \tab evaluate a function of several images \cr
	\code{\link{scaletointerval}} \tab rescale pixel values \cr
	\code{\link{zapsmall.im}} \tab set very small pixel values to zero \cr
	\code{\link{levelset}} \tab level set of an image\cr 
	\code{\link{solutionset}} \tab region where an expression is true \cr
	\code{\link{imcov}} \tab spatial covariance function of image \cr
	\code{\link{convolve.im}} \tab spatial convolution of images \cr
	\code{\link{pixelcentres}} \tab extract centres of pixels \cr
	\code{\link{transmat}} \tab convert matrix of pixel values \cr
	                       \tab to a different indexing convention
     }

    \bold{Line segment patterns}

    An object of class \code{"psp"} represents a pattern of straight line
    segments.
    \tabular{ll}{
      \code{\link{psp}} \tab create a line segment pattern \cr
      \code{\link{as.psp}} \tab convert other data into a line segment pattern \cr
      \code{\link{edges}} \tab extract edges of a window \cr
      \code{\link{is.psp}} \tab determine whether a dataset has class \code{"psp"} \cr
      \code{\link{plot.psp}} \tab plot a line segment pattern \cr
      \code{\link{print.psp}} \tab print basic information \cr
      \code{\link{summary.psp}} \tab print summary information \cr
      \code{\link{[.psp}} \tab extract a subset of a line segment pattern \cr
      \code{\link{subset.psp}} \tab extract subset of line segment pattern \cr
      \code{\link{as.data.frame.psp}} \tab
      convert line segment pattern to data frame \cr
      \code{\link{marks.psp}} \tab extract marks of line segments \cr
      \code{\link{marks<-.psp}} \tab assign new marks to line segments \cr
      \code{\link{unmark.psp}} \tab delete marks from line segments \cr
      \code{\link{midpoints.psp}} \tab
      compute the midpoints of line segments \cr
      \code{\link{endpoints.psp}} \tab extract the endpoints of line segments \cr
      \code{\link{lengths_psp}} \tab compute the lengths of line segments \cr
      \code{\link{angles.psp}} \tab compute the orientation angles of line segments \cr
      \code{\link{superimpose}} \tab combine several line segment patterns  \cr
      \code{\link{flipxy}} \tab swap \eqn{x} and \eqn{y} coordinates \cr      
      \code{\link{rotate.psp}} \tab rotate a line segment pattern \cr      
      \code{\link{shift.psp}} \tab shift a line segment pattern \cr      
      \code{\link{periodify}} \tab make several shifted copies \cr      
      \code{\link{affine.psp}} \tab apply an affine transformation \cr      
      \code{\link{pixellate.psp}} \tab approximate line segment pattern
    by pixel image \cr      
      \code{\link{as.mask.psp}} \tab approximate line segment pattern
    by binary mask \cr      
      \code{\link{distmap.psp}} \tab compute the distance map of a line
      segment pattern \cr
      \code{\link{distfun.psp}} \tab compute the distance map of a line
      segment pattern \cr
      \code{\link{selfcrossing.psp}} \tab find crossing points between
      line segments \cr
      \code{\link{selfcut.psp}} \tab cut segments where they cross \cr
      \code{\link{crossing.psp}} \tab find crossing points between
      two line segment patterns \cr
      \code{\link{extrapolate.psp}} \tab extrapolate line segments to
      infinite lines \cr
      \code{\link{nncross}} \tab find distance to nearest line segment
      from a given point\cr
      \code{\link{nearestsegment}} \tab find line segment closest to a
      given point \cr
      \code{\link{project2segment}} \tab find location along a line segment
      closest to a given point \cr
      \code{\link{pointsOnLines}} \tab generate points evenly spaced
      along line segment \cr
      \code{\link{rlinegrid}} \tab generate a random array of parallel
      lines through a window
    }

    \bold{Tessellations}

    An object of class \code{"tess"} represents a tessellation.

    \tabular{ll}{
      \code{\link{tess}} \tab create a tessellation \cr
      \code{\link{quadrats}} \tab create a tessellation of rectangles\cr
      \code{\link{hextess}} \tab create a tessellation of hexagons \cr
      \code{\link{polartess}} \tab tessellation using polar coordinates \cr
      \code{\link{quantess}} \tab quantile tessellation \cr
      \code{\link{venn.tess}} \tab Venn diagram tessellation \cr 
      \code{\link{dirichlet}} \tab compute Dirichlet-Voronoi tessellation of points\cr
      \code{\link{delaunay}} \tab compute Delaunay triangulation of points\cr
      \code{\link{as.tess}} \tab convert other data to a tessellation \cr
      \code{\link{plot.tess}} \tab plot a tessellation \cr
      \code{\link{tiles}} \tab extract all the tiles of a tessellation \cr
      \code{\link{[.tess}} \tab extract some tiles of a tessellation \cr
      \code{\link{[<-.tess}} \tab change some tiles of a tessellation \cr
      \code{\link{intersect.tess}} \tab intersect two tessellations \cr
                              \tab or restrict a tessellation to a window \cr
      \code{\link{chop.tess}} \tab subdivide a tessellation by a line \cr
      \code{\link{tile.areas}}	\tab
      area of each tile in tessellation \cr
      \code{\link{bdist.tiles}}	\tab
      boundary distance for each tile in tessellation \cr
      \code{\link{connected.tess}} \tab find connected components of tiles \cr
      \code{\link{shift.tess}} \tab shift a tessellation \cr
      \code{\link{rotate.tess}} \tab rotate a tessellation \cr
      \code{\link{reflect.tess}} \tab reflect about the origin \cr
      \code{\link{flipxy.tess}} \tab reflect about the diagonal \cr
      \code{\link{affine.tess}} \tab apply affine transformation 
    }

    \bold{Three-dimensional point patterns}

    An object of class \code{"pp3"} represents a three-dimensional
    point pattern in a rectangular box. The box is represented by
    an object of class \code{"box3"}.

    \tabular{ll}{
      \code{\link{pp3}} \tab create a 3-D point pattern \cr
      \code{\link{plot.pp3}} \tab plot a 3-D point pattern \cr
      \code{\link{coords}} \tab extract coordinates \cr
      \code{\link{as.hyperframe}} \tab extract coordinates \cr
      \code{\link{subset.pp3}} \tab extract subset of 3-D point pattern \cr
      \code{\link{unitname.pp3}} \tab name of unit of length \cr
      \code{\link{npoints}} \tab count the number of points  \cr
      \code{\link{box3}} \tab create a 3-D rectangular box \cr
      \code{\link{as.box3}} \tab convert data to 3-D rectangular box \cr
      \code{\link{unitname.box3}} \tab name of unit of length \cr
      \code{\link{diameter.box3}} \tab diameter of box \cr
      \code{\link{volume.box3}} \tab volume of box \cr
      \code{\link{shortside.box3}} \tab shortest side of box \cr
      \code{\link{eroded.volumes}} \tab volumes of erosions of box 
    }

    \bold{Multi-dimensional space-time point patterns}

    An object of class \code{"ppx"} represents a 
    point pattern in multi-dimensional space and/or time.

    \tabular{ll}{
      \code{\link{ppx}} \tab create a multidimensional space-time point pattern \cr
      \code{\link{coords}} \tab extract coordinates \cr
      \code{\link{as.hyperframe}} \tab extract coordinates \cr
      \code{\link{subset.ppx}} \tab extract subset \cr
      \code{\link{unitname.ppx}} \tab name of unit of length \cr
      \code{\link{npoints}} \tab count the number of points  \cr
      \code{\link{boxx}} \tab define multidimensional box  \cr
      \code{\link{diameter.boxx}} \tab diameter of box \cr
      \code{\link{volume.boxx}} \tab volume of box \cr
      \code{\link{shortside.boxx}} \tab shortest side of box \cr
      \code{\link{eroded.volumes.boxx}} \tab volumes of erosions of box 
    }
    
    \bold{Linear networks}

    An object of class \code{"linnet"} represents a linear network
    (for example, a road network). This is supported in the
    sub-package \pkg{spatstat.linnet}.

    An object of class \code{"lpp"} represents a 
    point pattern on a linear network (for example,
    road accidents on a road network).
    
    \bold{Hyperframes}

    A hyperframe is like a data frame, except that the entries
    may be objects of any kind.

    \tabular{ll}{
      \code{\link{hyperframe}} \tab create a hyperframe \cr
      \code{\link{as.hyperframe}} \tab convert data to hyperframe \cr
      \code{\link{plot.hyperframe}} \tab plot hyperframe \cr
      \code{\link{with.hyperframe}} \tab evaluate expression using each row
      of hyperframe \cr
      \code{\link{cbind.hyperframe}} \tab combine hyperframes by columns\cr
      \code{\link{rbind.hyperframe}} \tab combine hyperframes by rows\cr
      \code{\link{as.data.frame.hyperframe}} \tab convert hyperframe to
      data frame \cr
      \code{\link{subset.hyperframe}} \tab method for \code{subset} \cr
      \code{\link{head.hyperframe}} \tab first few rows of hyperframe \cr
      \code{\link{tail.hyperframe}} \tab last few rows of hyperframe
    }
    
    \bold{Layered objects}

    A layered object represents data that should be plotted in
    successive layers, for example, a background and a foreground.

     \tabular{ll}{
      \code{\link{layered}} \tab create layered object \cr
      \code{\link{plot.layered}} \tab plot layered object\cr
      \code{\link{[.layered}} \tab extract subset of layered object
    }
    
    \bold{Colour maps}

    A colour map is a mechanism for associating colours with data.
    It can be regarded as a function, mapping data to colours.
    Using a \code{colourmap} object in a plot command
    ensures that the mapping from numbers to colours is
    the same in different plots. 

     \tabular{ll}{
      \code{\link{colourmap}} \tab create a colour map  \cr
      \code{\link{plot.colourmap}} \tab plot the colour map only\cr
      \code{\link{tweak.colourmap}} \tab alter individual colour values \cr
      \code{\link{interp.colourmap}} \tab make a smooth transition
      between colours \cr
      \code{\link{beachcolourmap}} \tab one special colour map
    }
  
  \bold{Inspection of data:}
  \tabular{ll}{
    \code{\link{summary}(X)} \tab
    print useful summary of point pattern \code{X}\cr
    \code{X} \tab
    print basic description of point pattern \code{X}  \cr
    \code{any(duplicated(X))} \tab
    check for duplicated points in pattern \code{X} \cr
    \code{\link{intensity}} \tab Mean intensity \cr
    \code{\link{quadratcount}} \tab Quadrat counts \cr
  }

  \bold{Distances in a point pattern:}

  \tabular{ll}{
    \code{\link{nndist}} \tab nearest neighbour distances \cr
    \code{\link{nnwhich}} \tab find nearest neighbours \cr
    \code{\link{pairdist}} \tab distances between all pairs of points\cr
    \code{\link{crossdist}} \tab distances between points in two patterns\cr
    \code{\link{nncross}} \tab nearest neighbours between two point patterns \cr
    \code{\link{exactdt}} \tab distance from any location to nearest data point\cr 
    \code{\link{distmap}} \tab distance map image\cr
    \code{\link{distfun}} \tab distance map function\cr
    \code{\link{nnmap}} \tab nearest point image \cr
    \code{\link{nnfun}} \tab nearest point function 
 }

  \bold{Programming tools:}
  \tabular{ll}{
    \code{\link{applynbd}} \tab apply function to every neighbourhood
    in a point pattern \cr
    \code{\link{markstat}} \tab apply function to the marks of neighbours
    in a point pattern \cr
    \code{\link{pppdist}} \tab find the optimal match between two point
    patterns
  }

  \bold{Distances in a three-dimensional point pattern:}

  \tabular{ll}{
    \code{\link{pairdist.pp3}} \tab distances between all pairs of
    points \cr
    \code{\link{crossdist.pp3}} \tab distances between points in
    two patterns \cr
    \code{\link{nndist.pp3}} \tab nearest neighbour distances \cr
    \code{\link{nnwhich.pp3}} \tab find nearest neighbours \cr
    \code{\link{nncross.pp3}} \tab find nearest neighbours in another pattern
  }

  \bold{Distances in multi-dimensional point pattern:}

  These are for multi-dimensional space-time
  point pattern objects (class \code{ppx}).

  \tabular{ll}{
    \code{\link{pairdist.ppx}} \tab distances between all pairs of
    points \cr
    \code{\link{crossdist.ppx}} \tab distances between points in
    two patterns \cr
    \code{\link{nndist.ppx}} \tab nearest neighbour distances \cr
    \code{\link{nnwhich.ppx}} \tab find nearest neighbours
  }

}

\section{Licence}{
  This library and its documentation are usable under the terms of the "GNU 
  General Public License", a copy of which is distributed with the package.
}
\author{
  \spatstatAuthors.
}
\section{Acknowledgements}{
  Kasper Klitgaard Berthelsen,
  Ottmar Cronie,
  Tilman Davies,
  Yongtao Guan,
  Ute Hahn,
  Abdollah Jalilian,
  Marie-Colette van Lieshout,
  Greg McSwiggan,
  Tuomas Rajala,
  Suman Rakshit,
  Dominic Schuhmacher,
  Rasmus Waagepetersen
  and
  Hangsheng Wang
  made substantial contributions of code.

  Additional contributions and suggestions from
  Monsuru Adepeju,
  Corey Anderson,
  Ang Qi Wei,
  Ryan Arellano,
  Jens \ifelse{latex}{\out{{\AA}str{\" o}m}}{Astrom},
  Robert Aue,
  Marcel Austenfeld,
  Sandro Azaele,
  Malissa Baddeley,
  Guy Bayegnak,
  Colin Beale,
  Melanie Bell,
  Thomas Bendtsen,
  Ricardo Bernhardt,
  Andrew Bevan,
  Brad Biggerstaff,
  Anders Bilgrau,
  Leanne Bischof,
  Christophe Biscio,
  Roger Bivand,
  Jose M. Blanco Moreno,
  Florent Bonneu,
  Jordan Brown,
  Ian Buller,
  Julian Burgos,
  Simon Byers,
  Ya-Mei Chang,
  Jianbao Chen,
  Igor Chernayavsky,
  Y.C. Chin,
  Bjarke Christensen,
  \ifelse{latex}{\out{Luc\'{\i{}}a}}{Lucia} Cobo Sanchez,
  Jean-Francois Coeurjolly,
  Kim Colyvas,
  Hadrien Commenges,
  Rochelle Constantine,
  Robin Corria Ainslie,
  Richard Cotton,
  Marcelino de la Cruz,
  Peter Dalgaard,
  Mario D'Antuono,
  Sourav Das,
  Peter Diggle,
  Patrick Donnelly,
  Ian Dryden,
  Stephen Eglen,
  Ahmed El-Gabbas,
  Belarmain Fandohan,
  Olivier Flores,
  David Ford,
  Peter Forbes,
  Shane Frank,
  Janet Franklin,
  Funwi-Gabga Neba,
  Oscar Garcia,
  Agnes Gault,
  Jonas Geldmann,
  Marc Genton,
  Shaaban Ghalandarayeshi,
  Julian Gilbey,
  Jason Goldstick,
  Pavel Grabarnik,
  C. Graf,
  Ute Hahn,
  Andrew Hardegen,
  Martin \Bogsted Hansen,
  Martin Hazelton,
  Juha Heikkinen,
  Mandy Hering,
  Markus Herrmann,
  Maximilian Hesselbarth,
  Paul Hewson,
  Hamidreza Heydarian,
  Kassel Hingee,
  Kurt Hornik,
  Philipp Hunziker,
  Jack Hywood,
  Ross Ihaka,
  \ifelse{latex}{\out{\u{C}enk I\c{c}\"{o}s}}{Cenk Icos},
  Aruna Jammalamadaka,
  Robert John-Chandran,
  Devin Johnson,
  Mahdieh Khanmohammadi,
  Bob Klaver,
  Lily Kozmian-Ledward,
  Peter Kovesi,
  Mike Kuhn,
  Jeff Laake,
  Robert Lamb,
  \ifelse{latex}{\out{Fr\'{e}d\'{e}ric}}{Frederic} Lavancier,
  Tom Lawrence,
  Tomas Lazauskas,
  Jonathan Lee,
  George Leser,
  Angela Li,
  Li Haitao,
  George Limitsios,
  Andrew Lister,
  Nestor Luambua,
  Ben Madin,
  Martin Maechler,
  Kiran Marchikanti,
  Jeff Marcus,
  Robert Mark,
  Peter McCullagh,
  Monia Mahling,
  Jorge Mateu Mahiques,
  Ulf Mehlig,
  Frederico Mestre,
  Sebastian Wastl Meyer,
  Mi Xiangcheng,
  Lore De Middeleer,
  Robin Milne,
  Enrique Miranda,
  Jesper \Moller,
  Annie \ifelse{latex}{\out{Molli{\'e}}}{Mollie},
  Ines Moncada,
  Mehdi Moradi,
  Virginia Morera Pujol,
  Erika Mudrak,
  Gopalan Nair,
  Nader Najari,
  Nicoletta Nava,
  Linda Stougaard Nielsen,
  Felipe Nunes,
  Jens Randel Nyengaard,
  Jens \Oehlschlaegel,
  Thierry Onkelinx,
  Sean O'Riordan,
  Evgeni Parilov,
  Jeff Picka,
  Nicolas Picard,
  Tim Pollington,
  Mike Porter,
  Sergiy Protsiv,
  Adrian Raftery,
  Suman Rakshit,
  Ben Ramage,
  Pablo Ramon,
  Xavier Raynaud,
  Nicholas Read,
  Matt Reiter,
  Ian Renner,
  Tom Richardson,
  Brian Ripley,
  Ted Rosenbaum,
  Barry Rowlingson,
  Jason Rudokas,
  Tyler Rudolph,
  John Rudge,
  Christopher Ryan,
  Farzaneh Safavimanesh,
  Aila \Sarkka,
  Cody Schank,
  Katja Schladitz,
  Sebastian Schutte,
  Bryan Scott,
  Olivia Semboli,
  \ifelse{latex}{\out{Fran\c{c}ois S\'{e}m\'{e}curbe}}{Francois Semecurbe},
  Vadim Shcherbakov,
  Shen Guochun,
  Shi Peijian,
  Harold-Jeffrey Ship,
  Tammy L Silva,
  Ida-Maria Sintorn,
  Yong Song, 
  Malte Spiess,
  Mark Stevenson,
  Kaspar Stucki,
  Jan Sulavik,
  Michael Sumner,
  P. Surovy,
  Ben Taylor,
  Thordis Linda Thorarinsdottir,
  Leigh Torres,
  Berwin Turlach,
  Torben Tvedebrink,
  Kevin Ummer,
  Medha Uppala,
  Andrew van Burgel,
  Tobias Verbeke,
  Mikko Vihtakari,
  Alexendre Villers,
  Fabrice Vinatier,
  Maximilian Vogtland,
  Sasha Voss,
  Sven Wagner,
  Hao Wang,
  H. Wendrock,
  Jan Wild,
  Carl G. Witthoft,
  Selene Wong,
  Maxime Woringer,
  Luke Yates,
  Mike Zamboni
  and
  Achim Zeileis.
}
\keyword{spatial}
\keyword{package}

