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

  An object of class \code{"metric"} defines a metric,
  and also supports many computations that involve distances measured
  using the metric.

  An object of class \code{"metric"} is essentially a list of functions,
  including at least the following elements:

  \tabular{ll}{
    \code{pairdist} \tab analogue of \code{\link{pairdist.ppp}} \cr
    \code{crossdist} \tab analogue of \code{\link{crossdist.ppp}} \cr
    \code{nndist} \tab analogue of \code{\link{nndist.ppp}} \cr
    \code{nnwhich} \tab analogue of \code{\link{nnwhich.ppp}} \cr
    \code{nncross} \tab analogue of \code{\link{nncross.ppp}} \cr
    \code{distmap} \tab analogue of \code{\link{distmap.owin}} 
  }

  Objects of class \code{"metric"} are produced by the function
  \code{\link{convexmetric}}.
}
\seealso{
  \code{\link{convexmetric}}
}
\examples{
  K <- owin(poly=list(x=c(2.5,2,0.5,-2.5,-2,-0.5),y=c(0,1,2,0,-1,-2)))
  m <- convexmetric(K)
  m
  (m[["nndist"]])(cells)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{attribute}
 
 
