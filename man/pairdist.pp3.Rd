\name{pairdist.pp3}
\alias{pairdist.pp3}
\title{Pairwise distances in Three Dimensions}
\description{
  Computes the matrix of distances between all pairs of points
  in a three-dimensional point pattern.
}
\usage{
  \method{pairdist}{pp3}(X, \dots, periodic=FALSE, squared=FALSE)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"pp3"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{periodic}{
    Logical. Specifies whether to apply a periodic edge correction.
  }
  \item{squared}{
    Logical. If \code{squared=TRUE}, the squared distances are
    returned instead (this computation is faster).
  }
}
\value{
  A square matrix whose \code{[i,j]} entry is the distance
  between the points numbered \code{i} and \code{j}.
}
\details{
  This is a method for the generic function \code{pairdist}.

  Given a three-dimensional point pattern \code{X}
  (an object of class \code{"pp3"}),
  this function computes the Euclidean distances between all pairs of
  points in \code{X}, and returns the matrix of distances.

  Alternatively if \code{periodic=TRUE} and the window containing \code{X} is a
  box, then the distances will be computed in the \sQuote{periodic}
  sense (also known as \sQuote{torus} distance): opposite faces of the
  box are regarded as equivalent.
  This is meaningless if the window is not a box.

  If \code{squared=TRUE} then the \emph{squared} Euclidean distances
  \eqn{d^2} are returned, instead of the Euclidean distances \eqn{d}.
  The squared distances are faster to calculate, and are sufficient for
  many purposes (such as finding the nearest neighbour of a point).
}
\section{Distance values}{
  The values returned by \code{pairdist(X)} are distances,
  expressed as multiples of the unit of length of the spatial
  coordinates in \code{X}.
  The unit of length is given by \code{\link{unitname}(X)}.

  Note that, if the unit of length in \code{X} is a composite
  expression such as \sQuote{2 microns},
  then the values of \code{pairdist(X)} are expressed as multiples
  of 2 microns, rather than being expressed in microns.
}
\seealso{
  \code{\link{pairdist}},
  \code{\link{crossdist}},
  \code{\link{nndist}},
  \code{\link[spatstat.explore]{K3est}}
}
\examples{
if(require(spatstat.random)) {
   X <- runifpoint3(20)
} else {
   X <- osteo$pts[[1]]
}
   d <- pairdist(X)
   d <- pairdist(X, periodic=TRUE)
   d <- pairdist(X, squared=TRUE)
}
\author{
  \adrian
  based on two-dimensional code by 
  \pavel.
}
\keyword{spatial}
\keyword{math}

\concept{Three-dimensional}
