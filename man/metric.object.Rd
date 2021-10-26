\name{metric.object}
\alias{metric.object} %DoNotExport
\title{Distance Metric}
\description{
  An object of class \code{"metric"} defines a measure of distance
  between points, and supports many operations that involve distances.
}
\details{
  A \sQuote{metric} \eqn{d} is a measure of distance between points that
  satisfies
  \enumerate{
    \item
    \eqn{d(x,x) = 0} for any point \eqn{x},
    \item
    \eqn{d(x,y) > 0} for any two distinct points \eqn{x} and \eqn{y}
    \item
    symmetry: \eqn{d(x,y) = d(y,x)} for any two points \eqn{x} and \eqn{y}
    \item
    triangle inequality:
    \eqn{d(x,y) \le d(x,z) + d(z,y)}{d(x,y) <= d(x,z) + d(z,y)}
    for any three points \eqn{x,y,z}.
  }
  The Euclidean distance between points is an example of a metric.

  An object of class \code{"metric"} is a structure
  that defines a metric and supports many computations that involve
  the metric. The internal structure of this object,
  and the mechanism for performing these computations,
  are under development.

  Objects of class \code{"metric"} are produced by the function
  \code{\link{convexmetric}} and possibly by other functions.

  There are methods for \code{print} and \code{summary} 
  for the class \code{"metric"}.
  The \code{summary} method lists the operations that are supported
  by the metric.

  To perform distance calculations (for example, nearest-neighbour
  distances) using a desired metric instead of the
  Euclidean metric, first check whether the standard function for this
  purpose (for example \code{nndist.ppp}) 
  has an argument named \code{metric}. 
  the low-level function \code{\link{invoke.metric}}.
  If so, use the standard function
  and add the argument \code{metric}; if not, use
  the low-level function \code{\link{invoke.metric}}.
}
\seealso{
  \code{\link{convexmetric}}, \code{\link{invoke.metric}}
}
\examples{
  m <- convexmetric(square(c(-1,1)))
  summary(m)
  y <- nndist(cells, metric=m)
}

\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{attribute}
 
 
