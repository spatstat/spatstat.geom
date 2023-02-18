\name{integral.density}
\alias{integral.density}
\title{
  Compute Integral of One-Dimensional Kernel Density Estimate.
}
\description{
  Compute the integral of a kernel density estimate
  over a specified range.
}
\usage{
\method{integral}{density}(f, domain = NULL, weight=NULL, ...)
}
\arguments{
  \item{f}{
    A one-dimensional probability density estimate
    (object of class \code{"density"}) obtained from the
    function \code{\link[stats]{density.default}} or
    from \code{\link{unnormdensity}}.
  }
  \item{domain}{
    Optional. Range of values of the argument \eqn{x} over which the
    density \eqn{f(x)} should be integrated. A numeric vector of length
    2 giving the minimum and maximum values of \eqn{x}. Infinite limits
    are permitted.
  }
  \item{weight}{
    Optional. A \code{function(x)} specifying a weight integrand.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  This is a method for the generic function \code{\link{integral}}.
  It computes the numerical integral
  \deqn{
    I = \int f(x) dx
  }{
    I = integral f(x) dx
  }
  of the density estimate \code{f}.
  If \code{weight} is specified, then the weighted integral
  \deqn{
    I = \int w(x) f(x) dx
  }{
    I = integral w(x) f(x) dx
  }
  is computed, where \eqn{w} is the function specified by \code{weight}.
  This function must return finite numerical values.

  If \code{domain} is specified, the integral is restricted to the
  interval of \eqn{x} values given by the \code{domain}.

  Integrals are calculated numerically using the trapezoidal rule
  restricted to the domain given.
}
\value{
  A single numerical value.
}
\author{
  \adrian.
}
\seealso{
  \code{\link[stats]{density.default}}
  
  \code{\link[spatstat.explore]{quantile.density}}, \code{\link[spatstat.explore]{CDF.density}}
}
\examples{
  x <- runif(10)
  d <- density(x, bw=0.1)
  integral(d) # should be approximately 1
  integral(d, domain=c(-Inf, 0)) # mass on negative half-line
  ## mean of density
  integral(d, weight=function(x) x)
}
\keyword{methods}
\keyword{univar}
\keyword{nonparametric}
\keyword{math}
