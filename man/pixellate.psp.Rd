\name{pixellate.psp}  
\alias{pixellate.psp}
\title{
  Convert Line Segment Pattern to Pixel Image
}
\description{
  Converts a line segment pattern to a pixel image
  by measuring the length or number of lines intersecting each pixel.
}
\usage{
\method{pixellate}{psp}(x, W=NULL, ..., weights = NULL,
                        what=c("length", "number", "indicator"),
                        DivideByPixelArea=FALSE)
}
\arguments{
  \item{x}{
    Line segment pattern (object of class \code{"psp"}).
  }
  \item{W}{
    Optional window (object of class \code{"owin"}) determining
    the pixel resolution.
  }
  \item{\dots}{
    Optional arguments passed to \code{\link[spatstat.geom]{as.mask}} to determine
    the pixel resolution.
  }
  \item{weights}{
    Optional vector of weights associated with each line segment,
    used when \code{what="length"}.
  }
  \item{what}{
    String (partially matched) indicating whether to compute the
    total length of intersection of segments with each pixel
    (\code{what="length"}, the default),
    the total number of segments intersecting each pixel
    (\code{what="number"}) or an indicator variable that equals 1 if
    any segments intersect the pixel, and 0 otherwise
    (\code{code="indicator"}).
  }
  \item{DivideByPixelArea}{
    Logical value, indicating whether the resulting pixel values
    should be divided by the pixel area.
  }
}
\details{
  This function converts a line segment pattern to a pixel image by computing,
  for each pixel, the total length of
  intersection between the pixel and the line segments.
  Alternatively it can count the number of line segments intersecting
  each pixel. 
  
  This is a method for the generic function
  \code{\link{pixellate}} for the class of line segment patterns.

  The pixel raster is determined by \code{W}
  and the optional arguments \code{\dots}.
  If \code{W} is missing or \code{NULL}, it defaults to the window
  containing \code{x}.
  Then \code{W} is converted to a
  binary pixel mask using \code{\link[spatstat.geom]{as.mask}}. The arguments
  \code{\dots} are passed to \code{\link[spatstat.geom]{as.mask}} to
  control the pixel resolution.

  If \code{weights} are given and \code{what="length"},
  then the length of the intersection
  between line segment \code{i} and pixel \code{j} is multiplied by
  \code{weights[i]} before the lengths are summed for each pixel.
}
\value{
  A pixel image (object of class \code{"im"}) with numeric values.
}
\seealso{
  \code{\link{pixellate}},
  \code{\link[spatstat.geom]{as.mask}},
  \code{\link{psp2mask}}.
  
  Use \code{\link{psp2mask}} if you only want to know
  which pixels are intersected by lines.
}
\examples{
  X <- psp(runif(10),runif(10), runif(10), runif(10), window=owin())
  plot(pixellate(X))
  plot(X, add=TRUE)
  sum(lengths_psp(X))
  sum(pixellate(X))
  plot(pixellate(X, what="n"))
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{manip}

