\name{npoints}
\alias{npoints}
\alias{npoints.ppp}
\alias{npoints.pp3}
\alias{npoints.ppx}
\title{Number of Points in a Point Pattern}
\description{
  Returns the number of points in a point pattern of any kind.
}
\usage{
  npoints(x)
  \method{npoints}{ppp}(x)
  \method{npoints}{pp3}(x)
  \method{npoints}{ppx}(x)
}
\arguments{
  \item{x}{
    A point pattern (object of class \code{"ppp"},
    \code{"pp3"}, \code{"ppx"} or some other suitable class).
  }
}
\value{
  Integer.
}
\details{
  This function returns the number of points in
  a point pattern. The function \code{npoints} is generic
  with methods for the classes \code{"ppp"}, \code{"pp3"},
  \code{"ppx"} and possibly other classes.
}
\seealso{
  \code{\link{ppp.object}},
  \code{\link{print.pp3}},
  \code{\link{print.ppx}}.
}
\examples{
   npoints(cells)
}
\author{\adrian
  and \rolf
}
\keyword{spatial}
\keyword{manip}
