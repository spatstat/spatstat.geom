\name{as.function.tess}
\alias{as.function.tess}
\title{
  Convert a Tessellation to a Function
}
\description{
  Convert a tessellation into a function of the \eqn{x} and \eqn{y} coordinates.
  The default function values are factor levels specifying which tile of the
  tessellation contains the point \eqn{(x,y)}.
}
\usage{
  \method{as.function}{tess}(x,\dots,values=NULL)
}
\arguments{
  \item{x}{
    A tessellation (object of class \code{"tess"}).
  }
  \item{values}{
    Optional. A vector giving the values of the function for each tile
    of \code{x}.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  This command converts a tessellation
  (object of class \code{"tess"}) to a \code{function(x,y)}
  where the arguments \code{x} and \code{y} are (vectors of) spatial
  coordinates. The corresponding function values are factor levels
  identifying which tile of the tessellation contains each point.
  Values are \code{NA} if the corresponding point lies outside the
  tessellation.

  If the argument \code{values} is given, then it determines the value
  of the function in each tile of \code{x}.
}
\value{
  A function in the \R language, also belonging to the
  class \code{"funxy"} and \code{"tessfun"}.

  The class \code{"tessfun"}
  has methods for \code{plot}, \code{print}, \code{as.tess}
  and \code{integral}.
}
\author{
  \spatstatAuthors
}
\seealso{
  \code{\link{integral.tessfun}} for integration of the function.
  
  \code{\link{tileindex}} for the low-level calculation of tile index.

  \code{\link{cut.ppp}} and \code{\link{split.ppp}} to
  divide up the points of a point pattern according to
  a tessellation.
}
\examples{
  X <- runifrect(7)
  V <- dirichlet(X)
  f <- as.function(V)
  f(0.1, 0.4)
  plot(f)
}
\keyword{spatial}
\keyword{manip}
\concept{Tessellation}
