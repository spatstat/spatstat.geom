\name{psp2mask}
\alias{psp2mask}
\alias{as.mask.psp}
\title{
  Convert Line Segment Pattern to Binary Pixel Mask
}
\description{
  Converts a line segment pattern to a binary pixel mask
  by determining which pixels intersect the lines.
}
\usage{
psp2mask(x, W=NULL, ...)
as.mask.psp(x, W=NULL, ...)
}

\arguments{
  \item{x}{
    Line segment pattern (object of class \code{"psp"}).
  }
  \item{W}{
    Optional window (object of class \code{"owin"})
    determining the pixel raster.
  }
  \item{\dots}{
    Optional extra arguments passed to \code{\link{as.mask}}
    to determine the pixel resolution.
  }
}
\details{
  The functions \code{psp2mask} and \code{as.mask.psp} are currently
  identical. In future versions of the package, \code{as.mask.psp}
  will be deprecated, and then removed. 
  
  This function converts a line segment pattern to a binary
  pixel mask by determining which pixels intersect the lines.

  The pixel raster is determined by \code{W}
  and the optional arguments \code{\dots}.
  If \code{W} is missing or \code{NULL}, it defaults to the window
  containing \code{x}.
  Then \code{W} is converted to a
  binary pixel mask using \code{\link{as.mask}}. The arguments
  \code{\dots} are passed to \code{\link{as.mask}} to
  control the pixel resolution.
}
\value{
  A window (object of class \code{"owin"})
  which is a binary pixel mask (type \code{"mask"}).
}
\seealso{
  \code{\link{pixellate.psp}},
  \code{\link{as.mask}}.

  Use \code{\link{pixellate.psp}} if you want to measure the
  length of line in each pixel.
}
\examples{
  X <- psp(runif(10), runif(10), runif(10), runif(10), window=owin())
  plot(psp2mask(X))
  plot(X, add=TRUE, col="red")
}
\author{
  \adrian
  and \rolf
}
\keyword{spatial}
\keyword{manip}

