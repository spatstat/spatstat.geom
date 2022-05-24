\name{mean.ewcdf} %DontDeclareMethods
\alias{mean.ewcdf}
\alias{mean.ecdf}
\title{Mean of Empirical Cumulative Distribution Function}
\description{
  Calculates the mean of a (weighted or unweighted)
  empirical cumulative distribution function.
}
\usage{
\method{mean}{ecdf}(x, trim=0, \dots)

\method{mean}{ewcdf}(x, trim=0, \dots)
}
\arguments{
  \item{x}{
    An empirical cumulative distribution function
    (object of class \code{"ecdf"} created by \code{\link[stats]{ecdf}})
    or a weighted empirical cumulative distribution function
    (object of class \code{"ewcdf"} created by \code{\link{ewcdf}}).
  }
  \item{trim}{
    The fraction (0 to 0.5) of data values to be trimmed from
    each end of their range, before the mean is computed.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  These functions are methods for the generic
  \code{\link[base]{mean}} 
  for the classes \code{"ecdf"} and \code{"ewcdf"}.

  They calculate the mean of the probability distribution
  corresponding to the cumulative distribution function \code{x}.
  This is equivalent to calculating the (weighted or unweighted)
  mean of the original data values.

  For \emph{weighted} empirical cumulative distribution functions
  (class \code{"ewcdf"}) the weights will first be normalised so that they
  sum to 1. The result of \code{mean.ewcdf}
  is always an average or weighted average or the original data values.
  The argument \code{trim} is interpreted as a probability
  under this normalised distribution; the corresponding
  quantiles are computed, and data outside these quantiles is deleted
  before calculating the weighted mean.
}
\value{
  A single number.
}
\seealso{
  Generic \code{\link[base]{mean}} and
  \code{\link[stats]{weighted.mean}}.

  \code{\link[stats]{ecdf}}, \code{\link{ewcdf}}
  to create the cumulative distribution functions.

  \code{\link[spatstat.explore]{stieltjes}} for integration with respect to
  a cumulative distribution function.
}
\examples{
  x <- 1:5
  mean(x)
  mean(ecdf(x))
  w <- 1:5
  mean(ewcdf(x, w))
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{methods}
\keyword{univar}
