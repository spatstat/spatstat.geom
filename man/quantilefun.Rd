\name{quantilefun}
\alias{quantilefun}
\alias{quantilefun.im}
\alias{quantilefun.ecdf}
\alias{quantilefun.ewcdf}
\title{
  Quantile Function
}
\description{
  Return the inverse function of a cumulative distribution function.
}
\usage{
  quantilefun(x, \dots)

  \method{quantilefun}{im}(x, \dots, type=1)
 
  \method{quantilefun}{ecdf}(x, \dots, type=1)

  \method{quantilefun}{ewcdf}(x, \dots, type=1)
}
\arguments{
  \item{x}{
    Data for which the quantile function should be calculated.
    Either an object containing data (of class \code{"im"} or other
    classes) or an object representing a cumulative distribution function
    (of class \code{"ecdf"} or \code{"ewcdf"}).
  }
  \item{\dots}{
    Other arguments passed to methods.
  }
  \item{type}{
    Integer specifying the type of quantiles,
    as explained in \code{\link[stats]{quantile.default}}.
    Only types 1 and 2 are currently implemented.
  }
}
\details{
  Whereas the command \code{\link[stats]{quantile}} calculates
  the quantiles of a dataset corresponding to desired probabilities
  \eqn{p}, the command \code{quantilefun}
  returns a function which can be used to compute any quantiles of the
  dataset. 

  If \code{f <- quantilefun(x)} then \code{f} is a function such that
  \code{f(p)} is the quantile associated with any given probability \code{p}.
  For example \code{f(0.5)} is the median of the original data, and
  \code{f(0.99)} is the 99th percentile of the original data.

  If \code{x} is a pixel image (object of class \code{"im"})
  then the pixel values of \code{x} will be extracted
  and the quantile function of the pixel values is constructed.

  If \code{x} is an object representing a cumulative distribution
  function (object of class \code{"ecdf"} or \code{"ewcdf"}) then the
  quantile function of the original data is constructed.
}
\value{
  A function in the \R language.
}
\seealso{
  \code{\link{ewcdf}}, 
  \code{\link{quantile.ewcdf}}, 
  \code{\link[stats]{ecdf}},
  \code{\link[stats]{quantile}}
}
\examples{
  ## numeric data
  z <- rnorm(50)
  FZ <- ecdf(z)
  QZ <- quantilefun(FZ)
  QZ(0.5) # median value of z
  if(interactive()) plot(QZ, xlim=c(0,1), xlab="probability",
                     ylab="quantile of z")
  ## image data
  Z <- bei.extra$elev
  if(require(spatstat.explore)) {
    FE <- spatialcdf(Z, normalise=TRUE)
  } else {
    FE <- ecdf(Z[])
  }
  QE <- quantilefun(FE)
  QE(0.5) # median elevation
  if(interactive()) plot(QE, xlim=c(0,1),
     xlab="probability", ylab="quantile of elevation")
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{nonparametric}
