\name{distmap.psp}  
\alias{distmap.psp}
\title{
  Distance Map of Line Segment Pattern
}
\description{
  Computes the distance from each pixel to the nearest
  line segment in the given line segment pattern.
}
\usage{
  \method{distmap}{psp}(X, \dots, extras=TRUE, clip=FALSE, metric=NULL)
}
\arguments{
  \item{X}{A line segment pattern (object of class \code{"psp"}).
  }
  \item{\dots}{Arguments passed to \code{\link{as.mask}}
    to control pixel resolution.
  }
  \item{extras}{
    Logical value specifying whether to compute the
    additional attributes \code{"index"} and \code{"bdry"}
    described in Details.
  }
  \item{clip}{
    Logical value specifying whether the resulting pixel image
    should be clipped to the window of \code{X}.
  }
  \item{metric}{
    Optional. A distance metric
    (object of class \code{"metric"}, see \code{\link{metric.object}})
    which will be used to compute the distances.
  }
}
\value{
  A pixel image (object of class \code{"im"}) whose greyscale values
  are the values of the distance map.
  The return value has attributes \code{"index"} and \code{"bdry"}
  which are also pixel images. 
}
\details{
  The ``distance map'' of a line segment pattern \eqn{X} is the function
  \eqn{f} whose value \code{f(u)} is defined for any two-dimensional
  location \eqn{u} as the shortest distance from \eqn{u} to \eqn{X}.

  This function computes the distance map of the line segment pattern \code{X}
  and returns the distance map as a pixel image. The greyscale value
  at a pixel \eqn{u} equals the distance from \eqn{u}
  to the nearest line segment of the pattern \code{X}.
  Distances are computed using analytic geometry.

  The result is a pixel image. If \code{clip=TRUE} (the default), the
  pixel values are defined only inside the original window of \code{X},
  and are \code{NA} outside this window.
  If \code{clip=FALSE}, the pixel values are defined at every pixel in
  \code{Frame(X)}. Computation is faster when \code{clip=FALSE}.
  
  Additionally, if \code{extras=TRUE}, the return value 
  has two attributes, \code{"index"} and \code{"bdry"}, which are
  also pixel images. The grey values in \code{"bdry"} give the
  distance from each pixel to the bounding rectangle of the image.
  The grey values in \code{"index"} are integers identifying which
  line segment of \code{X} is closest.
  Computation is faster when \code{extras=FALSE}.
  
  This is a method for the generic function \code{\link{distmap}}.

  Note that this function gives the exact distance from the
  \emph{centre of each pixel} to the nearest line segment.
  To compute the exact distance from the points in a point pattern
  to the nearest line segment, use \code{\link{distfun}} or one of the
  low-level functions \code{\link{nncross}}
  or \code{\link{project2segment}}.
}
\seealso{
  \code{\link{distmap}},
  \code{\link{distmap.owin}},
  \code{\link{distmap.ppp}},
  \code{\link{distfun}},
  \code{\link{nncross}},
  \code{\link{nearestsegment}},
  \code{\link{project2segment}}.
}
\examples{
    a <- psp(runif(20),runif(20),runif(20),runif(20), window=owin())
    Z <- distmap(a)
    plot(Z)
    plot(a, add=TRUE)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{math}
