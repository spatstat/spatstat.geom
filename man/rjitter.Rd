\name{rjitter}
\alias{rjitter}
\alias{rjitter.ppp}
\title{Random Perturbation of a Point Pattern}
\description{
  Applies independent random displacements to each point in a point pattern.
}
\usage{

rjitter(X, \dots)

\method{rjitter}{ppp}(X, radius, retry=TRUE, giveup = 10000, trim=FALSE,
        \dots, nsim=1, drop=TRUE)
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
  }
  \item{retry}{
    What to do when a perturbed point lies outside the window
    of the original point pattern. If \code{retry=FALSE},
    the point will be lost; if \code{retry=TRUE},
    the algorithm will try again.
  }
  \item{giveup}{
    Maximum number of unsuccessful attempts.
  }
  \item{trim}{
    Logical value. If \code{TRUE}, the displacement radius
    for each data point will be constrained to be less than or equal to
    the distance from the data point to the window boundary.
    This ensures that all displaced points will fall inside the window.
  }
  \item{\dots}{Ignored.}
  \item{nsim}{Number of simulated realisations to be generated.}
  \item{drop}{
    Logical. If \code{nsim=1} and \code{drop=TRUE} (the default), the
    result will be a point pattern, rather than a list 
    containing a point pattern.
  }
}
\value{
  The result of \code{rjitter.ppp} is
  a point pattern (an object of class \code{"ppp"})
  or a list of point patterns.
}
\details{
  The function \code{rjitter} is generic, with methods for point
  patterns (described here) and for some other types of geometrical objects.

  Each of the points in the point pattern \code{X} is subjected to
  an independent random displacement. The displacement vectors are
  uniformly distributed in a circle of radius \code{radius}.

  If a displaced point lies outside the window, then if
  \code{retry=FALSE} the point will be lost. 

  However if \code{retry=TRUE}, the algorithm will try again: each time a
  perturbed point lies outside the window, the algorithm will reject it and
  generate another proposed perturbation of the original point,
  until one lies inside the window, or until \code{giveup} unsuccessful
  attempts have been made. In the latter case, any unresolved points
  will be included without any perturbation. The return value will
  always be a point pattern with the same number of points as \code{X}.

  If \code{trim=TRUE}, then the displacement radius for each data point
  will be constrained to be less than or equal to 
  the distance from the data point to the window boundary.
  This ensures that the randomly displaced points will
  always fall inside the window; no displaced points will be lost and no
  retrying will be required.
}
\examples{
   X <- rsyst(owin(), 10, 10)
   Y <- rjitter(X, 0.02)
   plot(Y)
   Z <- rjitter(X)
   U <- rjitter(X, 0.025, trim=TRUE)
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{rexplode}}
}
\keyword{spatial}
\keyword{datagen}
\keyword{manip}
