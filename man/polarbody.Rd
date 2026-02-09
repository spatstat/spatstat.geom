\name{polarbody}
\alias{polarbody}
\alias{voronoiFlower}
\title{
  Polar Body and Voronoi Flower
}
\description{
  Compute the polar body or the Voronoi flower of a two-dimensional domain.
}
\usage{
polarbody(W, origin = c(0, 0), npoly = 512)

voronoiFlower(W, origin = c(0, 0), npoly = 512)
}
\arguments{
  \item{W}{
    Window (object of class \code{"owin"}) representing the domain.
  }
  \item{origin}{
    The location of the origin for the support function.
    Either a numeric vector of length 2 giving the location,
    or a point pattern containing only one point,
    or a list with two entries named \code{x} and \code{y},
    or one of the character strings
    \code{"centroid"}, \code{"midpoint"},
    \code{"left"}, \code{"right"}, \code{"top"}, \code{"bottom"},
    \code{"topleft"}, \code{"bottomleft"}, \code{"topright"} or
    \code{"bottomright"} (partially matched).
  }
  \item{npoly}{
    Number of edges in the approximating polygon.
  }
}
\details{
  The polar body and Voronoi flower of a set \eqn{W}
  are geometrical constructions based on the support function
  \eqn{h_W(\theta)}{h[W](theta)} of \eqn{W}.
  
  The polar body of a two-dimensional domain \eqn{W}
  is the set whose boundary has polar coordinates
  \eqn{(\theta, 1/h_W(\theta))}{(theta, 1/h[W](theta))}
  for all directions \eqn{\theta}{theta}.

  The Voronoi flower is the set whose boundary has polar coordinates
  \eqn{(\theta, 2 h_W(\theta))}{(theta, 2 * h[W](theta))}.
}
\value{
  A polygonal window (object of class \code{"owin"})
  representing the polar body or the Voronoi flower.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{SupportFun}}
}
\examples{
  plot(polarbody(letterR, "midpoint"))

  plot(voronoiFlower(letterR, "midpoint"))
}
\keyword{spatial}
\keyword{manip}
\keyword{math}
