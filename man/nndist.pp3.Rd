\name{nndist.pp3}
\alias{nndist.pp3}
\title{Nearest neighbour distances in three dimensions}
\description{
  Computes the distance from each point to its nearest neighbour
  in a three-dimensional point pattern.
  Alternatively computes the distance to the
  second nearest neighbour, or third nearest, etc.
}
\usage{
  \method{nndist}{pp3}(X, \dots, k=1, by=NULL)
}
\arguments{
  \item{X}{
    Three-dimensional point pattern
    (object of class \code{"pp3"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{k}{
    Integer, or integer vector. The algorithm will compute the distance to the
    \code{k}th nearest neighbour. 
  }
  \item{by}{
    Optional. A factor, which separates \code{X} into groups.
    The algorithm will compute the distance to
    the nearest point in each group. 
  }
}
\value{
  Numeric vector or matrix containing the 
  nearest neighbour distances for each point.

  If \code{k = 1} (the default), the return value is a
  numeric vector \code{v} such that \code{v[i]} is the
  nearest neighbour distance for the \code{i}th data point.
  
  If \code{k} is a single integer, then the return value is a
  numeric vector \code{v} such that \code{v[i]} is the
  \code{k}th nearest neighbour distance for the
  \code{i}th data point.

  If \code{k} is a vector, then the return value is a
  matrix \code{m} such that \code{m[i,j]} is the
  \code{k[j]}th nearest neighbour distance for the
  \code{i}th data point.
}
\details{
  This function computes the Euclidean distance from each point
  in a three-dimensional
  point pattern to its nearest neighbour (the nearest other
  point of the pattern). If \code{k} is specified, it computes the
  distance to the \code{k}th nearest neighbour.

  The function \code{nndist} is generic; this function
  \code{nndist.pp3} is the method for the class \code{"pp3"}.

  The argument \code{k} may be a single integer, or an integer vector.
  If it is a vector, then the \eqn{k}th nearest neighbour distances are
  computed for each value of \eqn{k} specified in the vector.

  If there is only one point (if \code{x} has length 1),
  then a nearest neighbour distance of \code{Inf} is returned.
  If there are no points (if \code{x} has length zero)
  a numeric vector of length zero is returned.

  If the argument \code{by} is given, it should be a \code{factor},
  of length equal to the number of points in \code{X}.
  This factor effectively partitions \code{X} into subsets,
  each subset associated with one of the levels of \code{X}.
  The algorithm will then compute, for each point of \code{X},
  the distance to the nearest neighbour \emph{in each subset}.

  To identify \emph{which} point is the nearest neighbour of a given point,
  use \code{\link{nnwhich}}.

  To use the nearest neighbour distances for statistical inference,
  it is often advisable to use the edge-corrected empirical distribution,
  computed by \code{\link[spatstat.explore]{G3est}}.

  To find the nearest neighbour distances from one point pattern
  to another point pattern, use \code{\link{nncross}}.
}
\section{Distance values}{
  The values returned by \code{nndist(X)} are distances,
  expressed as multiples of the unit of length of the spatial
  coordinates in \code{X}.
  The unit of length is given by \code{\link{unitname}(X)}.

  Note that, if the unit of length in \code{X} is a composite
  expression such as \sQuote{2 microns},
  then the values of \code{nndist(X)} are expressed as multiples
  of 2 microns, rather than being expressed in microns.
}
\section{Warnings}{
  An infinite or \code{NA} value is returned if the
  distance is not defined (e.g. if there is only one point
  in the point pattern).
}
\seealso{
  \code{\link{nndist}},
  \code{\link{pairdist}},
  \code{\link[spatstat.explore]{G3est}},
  \code{\link{nnwhich}}
}
\examples{
   X <- pp3(runif(40), runif(40), runif(40), box3(c(0,1)))

   # nearest neighbours
   d <- nndist(X)

   # second nearest neighbours
   d2 <- nndist(X, k=2)

   # first, second and third nearest
   d1to3 <- nndist(X, k=1:3)

   # distance to nearest point in each group
   marks(X) <- factor(rep(letters[1:4], 10))
   dby <- nndist(X, by=marks(X))
}
\author{
  \adrian
  based on code for two dimensions by
  \pavel.
}
\keyword{spatial}
\keyword{math}

\concept{Three-dimensional}
