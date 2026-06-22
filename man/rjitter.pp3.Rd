\name{rjitter.pp3}
\alias{rjitter.pp3}
\title{Random Perturbation of a 3D Point Pattern}
\description{
  Applies independent random displacements to each point in a
  three-dimensional point pattern.
}
\usage{
\method{rjitter}{pp3}(X, radius, retry=TRUE, giveup = 10000, trim=FALSE,
        \dots, nsim=1, drop=TRUE, adjust=1)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
  }
  \item{radius}{
    Scale of perturbations. A positive numerical value.
    The displacement vectors will be uniformly
    distributed inside a sphere of this radius.
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
  \item{adjust}{
    Adjustment factor applied to the radius.
    A numeric value or numeric vector.
  }
}
\value{
  The result of \code{rjitter.pp3} is
  a three-dimensional point pattern (an object of class \code{"pp3"})
  or a list of 3D point patterns.

  Each point pattern has attributes \code{"radius"}
  and (if \code{retry=TRUE}) \code{"tries"}.
}
\details{
  The function \code{\link{rjitter}} is generic, with methods for point
  patterns and for some other types of geometrical objects.
  This help page describes the method for three-dimensional point
  patterns of class \code{"pp3"}.

  Each of the points in the point pattern \code{X} is subjected to
  an independent random displacement. The displacement vectors are
  uniformly distributed in a sphere of radius \code{radius}.

  If a displaced point lies outside the domain of the pattern, then if
  \code{retry=FALSE} the point will be lost. 

  However if \code{retry=TRUE}, the algorithm will try again: each time a
  perturbed point lies outside the domain, the algorithm will reject
  the perturbed point and
  generate another proposed perturbation of the original point,
  until one lies inside the domain, or until \code{giveup} unsuccessful
  attempts have been made. In the latter case, any unresolved points
  will be included, without any perturbation. The return value will
  always be a point pattern with the same number of points as \code{X}.

  If \code{trim=TRUE}, then the displacement radius for each data point
  will be constrained to be less than or equal to 
  the distance from the data point to the domain boundary.
  This ensures that the randomly displaced points will
  always fall inside the domain; no displaced points will be lost and no
  retrying will be required. However, it implies that a point lying
  exactly on the boundary will never be perturbed.

  If \code{adjust} is given, the jittering radius will be multiplied
  by \code{adjust}. This allows the user to specify
  that the radius should be a multiple of the default radius.

  The resulting point pattern
  has an attribute \code{"radius"} giving the value
  of \code{radius} used.
  If \code{retry=TRUE}, the resulting point pattern also has an attribute
  \code{"tries"} reporting the maximum number of trials needed to
  ensure that all jittered points were inside the domain.
}
\examples{
   X <- osteo$pts[[1]]
   Y <- rjitter(X, 5)
   plot(Y)
   Z <- rjitter(X)
   U <- rjitter(X, 5, trim=TRUE)
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{rjitter}}
}
\keyword{spatial}
\keyword{datagen}
\keyword{manip}
