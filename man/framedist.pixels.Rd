\name{framedist.pixels}
\alias{framedist.pixels}
\alias{framedist.pixels}
\title{Distance to Bounding Frame}
\description{
  Computes the distances from each pixel to the bounding rectangle.
}
\usage{
 framedist.pixels(w, \dots, style=c("image", "matrix", "coords"))
}
\arguments{
  \item{w}{A window (object of class \code{"owin"}).}
  \item{\dots}{
    Arguments passed to \code{\link[spatstat.geom]{as.mask}} to determine
    the pixel resolution.}
  \item{style}{
    Character string (partially matched) determining the format of
    the output: either \code{"matrix"}, \code{"coords"} or
    \code{"image"}.
  }
}
\value{
  If \code{style="image"}, a pixel image (object of class \code{"im"})
  containing the distances from each pixel in the image raster
  to the boundary of the window.

  If \code{style="matrix"},
  a matrix giving the distances from each pixel in the image raster
  to the boundary of the window. Rows of this matrix correspond to
  the \eqn{y} coordinate and columns to the \eqn{x} coordinate.
  
  If \code{style="coords"}, a list with three components
  \code{x,y,z}, where \code{x,y} are vectors of length \eqn{m,n}
  giving the \eqn{x} and \eqn{y} coordinates respectively,
  and \code{z} is an \eqn{m \times n}{m x n} matrix such that
  \code{z[i,j]} is the distance from \code{(x[i],y[j])} to the
  boundary of the window. Rows of this matrix correspond to the
  \eqn{x} coordinate and columns to the \eqn{y} coordinate.
  This result can be plotted with \code{persp}, \code{image}
  or \code{contour}.
}
\details{
  This function computes, for each pixel \eqn{u}
  in the rectangular frame \code{Frame(w)}, the shortest distance
  to the boundary of \eqn{Frame(w)}.

  The grid of pixels is determined by the arguments \code{"\dots"} 
  passed to \code{\link[spatstat.geom]{as.mask}}. The distance from each pixel to the
  boundary is calculated exactly, using analytic geometry.
}
\seealso{
  \code{\link{bdist.pixels}}.
}
\examples{
  opa <- par(mfrow=c(1,2))
  plot(framedist.pixels(letterR))
  plot(bdist.pixels(letterR))
  par(opa)
}
\author{
  \adrian
  and \rolf
}
\keyword{spatial}
\keyword{math}
