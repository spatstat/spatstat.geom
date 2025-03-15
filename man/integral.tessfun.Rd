\name{integral.tessfun}
\alias{integral.tessfun}
\title{
  Integrate a Function Which is Constant on Each Tile of a Tessellation
}
\description{
  Given a function which is constant on each tile of a tessellation,
  compute the integral of the function.
}
\usage{
 \method{integral}{tessfun}(f, domain = NULL, \dots)
}
\arguments{
  \item{f}{
    Integrand. A function of class \code{"tessfun"}
    (created by \code{\link[spatstat.geom]{as.function.tess}}).
  }
  \item{domain}{
    Optional window (object of class \code{"owin"})
    specifying a subregion to which the integral should be restricted.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  The command \code{\link[spatstat.univar]{integral}} is generic.
  This is the method for objects of class \code{"tessfun"}.
  
  The integrand \code{f} should be a function of class
  \code{"tessfun"} created by \code{\link[spatstat.geom]{as.function.tess}}.
  It represents a function which takes a constant value
  on each tile of a tessellation.

  The integral is calculated by multiplying the area of each tile
  by the value of the function in that tile, and summing.
  This avoids the need for discretisation and avoids concomitant
  discretisation errors.
}
\value{
  A single numeric value.
}
\author{
  \adrian.
}
\seealso{
  \code{\link[spatstat.geom]{as.function.tess}}.
}
\examples{
  V <- dirichlet(runifrect(5))
  f <- as.function(V, values=1/tile.areas(V))
  integral(f) # should be close to 5.
}
\keyword{spatial}
\keyword{math}
\concept{Tessellation}
