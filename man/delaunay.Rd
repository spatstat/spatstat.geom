\name{delaunay}
\alias{delaunay}
\title{Delaunay Triangulation of Point Pattern}
\description{
  Computes the Delaunay triangulation of a spatial point pattern.
}
\usage{
delaunay(X)
}
\arguments{
  \item{X}{Spatial point pattern (object of class \code{"ppp"}).}
}
\details{
  The Delaunay triangulation of a spatial point pattern \code{X}
  is defined as follows. First the Dirichlet/Voronoi tessellation
  based on \code{X} is computed; see \code{\link{dirichlet}}. This
  tessellation is extended to cover the entire two-dimensional plane.
  Then two points of \code{X}
  are defined to be Delaunay neighbours if their Dirichlet/Voronoi tiles
  share a common boundary. Every pair of Delaunay neighbours is
  joined by a straight line to make the Delaunay triangulation.
  The result is a tessellation, consisting of
  disjoint triangles. The union of these triangles is the convex hull of
  \code{X}.
}
\value{
  A tessellation (object of class \code{"tess"}). The window of the
  tessellation is the convex hull of \code{X}, not the original window
  of \code{X}.
}
\seealso{
  \code{\link{tess}},
  \code{\link{dirichlet}},
  \code{\link{convexhull.xy}},
  \code{\link{ppp}},
  \code{\link{delaunayDistance}},
  \code{delaunayNetwork}.
}
\examples{
  X <- runifrect(42)
  plot(delaunay(X))
  plot(X, add=TRUE)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{manip}
