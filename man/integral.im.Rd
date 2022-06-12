\name{integral.im}
\alias{integral}
\alias{integral.im}
\title{
  Integral of a Pixel Image
}
\description{
  Computes the integral of a pixel image.
}
\usage{
integral(f, domain=NULL, \dots)

\method{integral}{im}(f, domain=NULL, weight=NULL, \dots)
}
\arguments{
  \item{f}{
    A pixel image (object of class \code{"im"}) with pixel values
    that can be treated as numeric or complex values.
  }
  \item{domain}{
    Optional. Window specifying the domain of integration.
    Alternatively a tessellation.
  }
  \item{\dots}{
    Ignored.
  }
  \item{weight}{
   Optional. A pixel image (object of class \code{"im"})
   or a \code{function(x,y)} giving a numerical weight
   to be applied to the integration.
  }
 }
\details{
  The function \code{integral} is generic, with methods
  for \code{"im"}, \code{"msr"}, \code{"linim"} and \code{"linfun"}.
  
  The method \code{integral.im} treats the pixel image \code{f} as a function of
  the spatial coordinates, and computes its integral.
  The integral is calculated
  by summing the pixel values and multiplying by the area of one pixel.

  The pixel values of \code{f} may be numeric, integer, logical or
  complex. They cannot be factor or character values.

  The logical values \code{TRUE} and \code{FALSE} are converted to
  \code{1} and \code{0} respectively, so that the integral of a logical
  image is the total area of the \code{TRUE} pixels, in the same units
  as \code{unitname(x)}.

  If \code{domain} is a window (class \code{"owin"}) then the integration
  will be restricted to this window. If \code{domain} is a tessellation
  (class \code{"tess"}) then the integral of \code{f} in each
  tile of \code{domain} will be computed.

  If \code{weight} is given, it should be a pixel image or a function of
  coordinates \eqn{x} and \eqn{y} returning numerical values.
  Then each pixel value of \code{f}
  will be multiplied by the corresponding value of \code{weight}.
  Effectively, the result is the integral of \code{weight * f}.
}
\value{
  A single numeric or complex value (or a vector of such values
  if \code{domain} is a tessellation).
}
\seealso{
  \code{\link{eval.im}},
  \code{\link{[.im}}
}
\examples{
   # approximate integral of f(x,y) dx dy
   f <- function(x,y){3*x^2 + 2*y}
   Z <- as.im(f, square(1))
   integral(Z)
   # correct answer is 2

   # integrate over the subset [0.1,0.9] x [0.2,0.8]
   W <- owin(c(0.1,0.9), c(0.2,0.8))
   integral(Z, W)

   # weighted integral
   integral(Z, weight=function(x,y){x})
}
\author{
\spatstatAuthors.
}
\keyword{spatial}
\keyword{math}
