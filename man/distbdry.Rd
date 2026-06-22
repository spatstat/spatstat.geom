\name{distbdry}
\alias{distbdry}
\alias{distbdry.ppp}
\alias{distbdry.pp3}
\title{
  Distance To Boundary
}
\description{
  Compute the distance from each point in a point pattern
  to the boundary of the spatial domain.
}
\usage{
distbdry(X)

\method{distbdry}{ppp}(X)

\method{distbdry}{pp3}(X)
}
\arguments{
  \item{X}{
    A point pattern (class \code{"ppp"} or \code{"pp3"} or
    other classes).
  }
}
\details{
  For each point \code{X[i]} this function computes the shortest
  distance to the boundary of the observation window or domain of the
  point pattern.

  The function \code{distbdry} is generic, with methods for
  class \code{"ppp"}, \code{"pp3"} and other classes.
  
  For a two-dimensional spatial point pattern \code{X} of class \code{"ppp"},
  \code{distbdry(X)} is equivalent to \code{\link{bdist.points}(X)}.
}
\value{
  A numeric vector with one entry for each point of \code{X}.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{bdist.points}}
}
\examples{
 distbdry(cells)[1:10]
}
\keyword{spatial}
\keyword{math}
