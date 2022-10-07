\name{delaunayDistance}
\alias{delaunayDistance}
\title{Distance on Delaunay Triangulation}
\description{
  Computes the graph distance in the Delaunay triangulation
  of a point pattern.
}
\usage{
delaunayDistance(X)
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

  The \emph{graph distance} 
  in the Delaunay triangulation between two points \code{X[i]} and \code{X[j]}
  is the minimum number of edges of the Delaunay triangulation
  that must be traversed to go from \code{X[i]} to \code{X[j]}.
  Two points have graph distance 1 if they are immediate neighbours.

  This command returns a matrix \code{D} such that
  \code{D[i,j]} is the graph distance
  between \code{X[i]} and \code{X[j]}.
}
\section{Definition of neighbours}{
  Note that \code{\link{dirichlet}(X)}
  restricts the Dirichlet tessellation to the window containing
  \code{X}, whereas \code{dirichletDistance} uses the Dirichlet
  tessellation over the entire two-dimensional plane.
  Some points may be Delaunay neighbours
  according to \code{delaunayDistance(X)}
  although the corresponding tiles of \code{dirichlet(X)}
  do not share a boundary inside \code{Window(X)}.
}
\value{
  A symmetric square matrix with non-negative integer entries.
}
\seealso{
  \code{\link{delaunay}},
  \code{delaunayNetwork}.
}
\examples{
  X <- runifrect(20)
  M <- delaunayDistance(X)
  plot(delaunay(X), lty=3)
  text(X, labels=M[1, ], cex=2)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{manip}
