\name{simplify.owin}
\Rdversion{1.1}
\alias{simplify.owin}
\title{
  Approximate a Window by a Simpler Polygon
}
\description{
  Given a window, this function finds a simpler polygon
  that approximates it.
}
\usage{
simplify.owin(W, dmin)
}
\arguments{
  \item{W}{
    The window which is to be simplied.
    An object of class \code{"owin"}.
  }
  \item{dmin}{
    Numeric value. The smallest permissible length of an edge.
  }
}
\details{
  This function simplifies a polygon \code{W}
  by recursively deleting the shortest edge of \code{W}
  until all remaining edges are longer than the specified
  minimum length \code{dmin}, or until there are only three edges left.

  The argument \code{W} must be a window (object of class
  \code{"owin"}). If \code{W} is a rectangle, it is returned without alteration.
  If \code{W} is a binary mask, it is first converted to a polygonal
  window (with many small vertical and horizontal edges), then simplified.
}
\value{
  Another window (object of class \code{"owin"})
  of type \code{"polygonal"}.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{fillholes.owin}},
  \code{\link{owin}}
}
\examples{
  plot(letterR, col="red")
  plot(simplify.owin(letterR, 0.3), col="blue", add=TRUE)

  W <- Window(chorley)
  plot(W)
  WS <- simplify.owin(W, 2)
  plot(WS, add=TRUE, border="green")
  points(vertices(WS))
}
\keyword{spatial}
\keyword{math}
