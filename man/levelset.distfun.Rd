\name{levelset.distfun}
\alias{levelset.distfun}
\Rdversion{1.1}
\title{
  Level Set of a Distance Function
}
\description{
  Given a distance function, find the spatial region where
  the function value is less than a
  specified threshold value (or greater than a threshold, etc),
}
\usage{
\method{levelset}{distfun}(X, thresh, compare="<=", \dots)
}
\arguments{
  \item{X}{
    Object of class \code{"distfun"} representing the distance function
    of a spatial object.
  }
  \item{thresh}{
    Single numeric value of distance, defining the level set.
  }
  \item{compare}{Character string specifying one of the comparison
    operators \code{"<", ">", "==", "<=", ">=", "!="}. 
  }
  \item{\dots}{
    Arguments passed to \code{\link{as.mask}} if it is necessary to
    convert the result to a binary mask.
  }
}
\details{
  This is a method for the generic function \code{\link{levelset}}
  which computes level sets of a spatial variable.

  A distance function (object of class \code{"distfun"} created by
  \code{\link{distfun}}) is a function of spatial location
  which measures the distance to a given spatial object.

  If \code{f <- distfun(B)}, where \code{B} is some object,
  then the level set of \code{f} consisting of all spatial locations \code{u}
  where \code{f(u) <= thresh},
  is simply the dilation of the original object \code{B} by a distance equal
  to \code{thresh}. The command \code{levelset(X, thresh)} computes this
  dilation as a spatial window.

  In some cases the object must be converted to a binary pixel mask
  and the computation performed by \code{\link{levelset.im}}. In that
  case, any arguments \code{\dots} are passed to \code{\link{as.mask}}
  to determine the pixel resolution.
}
\value{
  A spatial window (object of class \code{"win"}).
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{distfun}},
  \code{\link{methods.distfun}},
  \code{\link{levelset}}
}
\examples{
  f <- distfun(cells)
  V <- levelset(f, 0.05)
  plot(V)
  plot(cells, add=TRUE)
}
\keyword{spatial}
\keyword{methods}
