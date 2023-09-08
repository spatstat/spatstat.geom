\name{integral}
\alias{integral}
\title{
  Integral of a Function or Spatial Object
}
\description{
  Computes the integral of a function or spatial object.
}
\usage{
integral(f, domain=NULL, \dots)
}
\arguments{
  \item{f}{
    A function, or a spatial object that can be treated as a function.
  }
  \item{domain}{
    Optional. Data specifying the domain of integration.
  }
  \item{\dots}{
    Arguments passed to methods.
  }
 }
\details{
  The function \code{integral} is generic.
  It calculates the integral of a function, or
  the integral of a spatial object that can be treated as a function.
  It has methods
  for one-dimensional functions (\code{"density"}, \code{"fv"})
  and for spatial objects
  (\code{"im"}, \code{"msr"}, \code{"linim"}, \code{"linfun"}).
}
\value{
  A single numeric or complex value, or a vector of such values.
}
\seealso{
  \code{\link[spatstat.geom]{integral.im}},
  \code{\link[spatstat.geom]{integral.density}}.
}
\author{
\spatstatAuthors.
}
\keyword{spatial}
\keyword{math}
