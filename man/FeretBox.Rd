\name{FeretBox}
\alias{FeretBox}
\title{
  Feret Box Enclosing a Domain
}
\description{
  Given a spatial domain, calculate its Feret box, the narrowest
  rectangle enclosing the domain at a given angle (or at any angle).
}
\usage{
FeretBox(W, theta = NULL)
}
\arguments{
  \item{W}{
    Window (object of class \code{"owin"}) representing the domain.
  }
  \item{theta}{
    Either a single numeric value specifying an angle
    such that the sides of the box should be parallel and perpendicular
    to \code{theta}, or \code{NULL} signifying that the orientation of
    the box is unconstrained.
  }
}
\details{
  If \code{theta} is given, then the smallest rectangle containing \code{W}
  with sides parallel and perpendicular to angle \code{theta}
  will be computed and returned.

  If \code{theta=NULL} (the default), then the orientation of the
  rectangle is unconstrained. The rectangle with the narrowest possible
  width (at any orientation) is determined by finding the minimum
  Feret diameter. 
}
\value{
  A polygonal window (object of class \code{"owin"}) representing a
  rotated rectangle.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{FeretDiamFun}}
}
\examples{
  W <- Window(chorley)
  B <- FeretBox(W)
  plot(B, main="", type="n")
  plot(W, add=TRUE)
  plot(B, add=TRUE, border="red")
}
\keyword{spatial}
\keyword{manip}
\keyword{math}
