\name{tile.centroids}
\alias{tile.centroids}
\title{Compute Centroids of Tiles in a Tessellation}
\description{
  Finds the centroid of each tile in a tessellation
  and returns them as a point pattern.
}
\usage{
tile.centroids(x)
}
\arguments{
  \item{x}{A tessellation (object of class \code{"tess"}).}
}
\details{
  A tessellation is a collection of disjoint spatial regions
  (called \emph{tiles}) that fit together to form a larger spatial
  region. See \code{\link{tess}}.

  This command determines the centroid of each of the tiles 
  that make up the tessellation \code{x} using \code{\link{centroid.owin}}.
  The result is a point pattern 
  containing the centroids,
  listed in the same order as the tiles would be listed by \code{tiles(x)}.
}
\value{
  A point pattern (object of class \code{"ppp"}).
}
\seealso{
  \code{\link{tess}},
  \code{\link{tile.areas}}
}
\examples{
  A <- tess(xgrid=0:2,ygrid=0:2)
  plot(A, main="")
  plot(tile.centroids(A), add=TRUE)
  D <- dirichlet(runifrect(6))
  plot(D)
  plot(tile.centroids(D), add=TRUE)
}
\author{\spatstatAuthors.}
\keyword{spatial}
\keyword{manip}
