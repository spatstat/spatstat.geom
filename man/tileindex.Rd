\name{tileindex}
\alias{tileindex}
\title{
  Determine Which Tile Contains Each Given Point
}
\description{
  Given a tessellation and a list of spatial points,
  determine which tile of the tessellation contains each of the
  given points.
}
\usage{
tileindex(x, y, Z, close.gaps=TRUE, all.inside=FALSE)
}
\arguments{
  \item{x,y}{
    Spatial coordinates.
    Numeric vectors of equal length.
    (Alternatively \code{y} may be missing
    and \code{x} may be an object containing
    spatial coordinates).
  }
  \item{Z}{
    A tessellation (object of class \code{"tess"}).
  }
  \item{close.gaps}{
    Logical value specifying whether all points \eqn{x,y}
    lying inside the window of \code{Z} must be classified
    as belonging to a tile. This avoids the effect of numerical errors
    in the tile geometry. See Details.
  }
  \item{all.inside}{
    Logical value specifying whether all points \eqn{x,y} should be 
    classified as lying inside a tile. This implies that points lying
    outside the window will be assigned to the closest tile. See Details.
  }
}
\details{
  This function determines which tile of the tessellation \code{Z}
  contains each of the spatial points
  with coordinates \code{(x[i],y[i])}.

  The result is a factor, of the same length as \code{x} and \code{y},
  indicating which tile contains each point. The levels of the factor
  are the names of the tiles of \code{Z}.

  A point lying outside the window containing the tessellation
  is assigned the value \code{NA} by default.
  However if \code{all.inside=TRUE} is specified,
  then every point will be assigned to
  the closest tile, and no \code{NA} values are returned.

  It is possible that, due to numerical error, a point lying \emph{inside} the
  window may not be classified as belonging to any of the tiles of
  \code{Z}. If this occurs, the default behaviour is to assign
  the point to the closest tile. This can be suppressed
  by setting \code{close.gaps=FALSE} (and \code{all.inside=FALSE});
  in that case, \code{NA} values are returned for such points.
}
\value{
  A factor, of the same length as \code{x} and \code{y},
  whose levels are the names of the tiles of \code{Z}.
  The factor values may include \code{NA} unless \code{all.inside=TRUE}
  was specified.
}
\author{
  \spatstatAuthors
}
\seealso{
  \code{\link{cut.ppp}} and \code{\link{split.ppp}} to
  divide up the points of a point pattern according to
  a tessellation.

  \code{\link{as.function.tess}} to create a function whose
  value is the tile index.
}
\examples{
  X <- runifrect(7)
  V <- dirichlet(X)
  tileindex(0.1, 0.4, V)
  tileindex(list(x=0.1, y=0.4), Z=V)
  tileindex(X, Z=V)
}
\keyword{spatial}
\keyword{manip}
