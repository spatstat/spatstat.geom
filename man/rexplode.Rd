\name{rexplode}
\alias{rexplode}
\alias{rexplode.ppp}
\title{
  Explode a Point Pattern by Displacing Duplicated Points
}
\description{
  Given a point pattern which contains duplicated points,
  separate the duplicated points from each other by
  slightly perturbing their positions.
}
\usage{
rexplode(X, \dots)

\method{rexplode}{ppp}(X, radius, \dots, nsim = 1, drop = TRUE)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
  }
  \item{radius}{
    Scale of perturbations. A positive numerical value.
    The displacement vectors will be uniformly
    distributed in a circle of this radius.
    There is a sensible default.
    Alternatively, \code{radius} may be a numeric vector of length
    equal to the number of points in \code{X}, giving a different
    displacement radius for each data point.
    Radii will be restricted to be less than or equal to the distance
    to the boundary of the window.
  }
  \item{\dots}{
    Ignored.
  }
  \item{nsim}{Number of simulated realisations to be generated.}
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE} (the default), the
    result will be a point pattern, rather than a list 
    containing a point pattern.
  }
}
\value{
  A point pattern (an object of class \code{"ppp"})
  or a list of point patterns.
}
\details{
  Duplicated points in the point pattern \code{X} are identified.
  Each group of duplicated points is then \sQuote{exploded} by
  randomly displacing the point locations to form a circular arrangement
  around the original position.
  
  This function is an alternative to \code{\link{rjitter.ppp}}.
  Whereas \code{\link{rjitter.ppp}}
  applies independent random displacements to each data point,
  \code{rexplode.ppp} applies displacements only to the points
  that are duplicated, and the displacements are mutually dependent
  within each group of duplicates, to ensure that the displaced points
  are visually separated from each other.

  First the code ensures that the displacement radius for each
  data point is less than or equal to the distance to the boundary
  of the window. Then each group of duplicated points (or data points
  with the same location but possibly different mark values) is
  taken in turn. The first element of the group is randomly displaced
  by a vector uniformly distributed in a circle of radius \code{radius}.
  The remaining elements of the group are then positioned around the
  original location, at the same distance from the orginal location,
  with equal angular spacing from the first point.
  The result is that each group of duplicated points becomes
  a circular pattern centred around the original location.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{rjitter.ppp}}
}
\examples{
   ## create a pattern containing duplicated points
   X <- runifrect(5) \%mark\% letters[1:5]
   X <- X[rep(1:5, 1 + rpois(5, 2))]
   ## explode it
   Y <- rexplode(X, 0.05)
   ## display
   if(interactive()) {
       plot(solist(X=X, 'explode(X)'=Y),
       main="", cols=2:6, cex=1.25, leg.side="bottom")
   }
}
\keyword{spatial}
\keyword{manip}
\keyword{datagen}
