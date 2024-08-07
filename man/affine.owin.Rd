\name{affine.owin} 
\alias{affine.owin} 
\title{Apply Affine Transformation To Window}
\description{
  Applies any affine transformation of the plane (linear transformation
  plus vector shift) to a window. 
}
\usage{
  \method{affine}{owin}(X, mat=diag(c(1,1)), vec=c(0,0), \dots, rescue=TRUE)
}
\arguments{
  \item{X}{Window (object of class \code{"owin"}).}
  \item{mat}{Matrix representing a linear transformation.}
  \item{vec}{Vector of length 2 representing a translation.}
  \item{rescue}{
    Logical. If \code{TRUE}, the transformed window
    will be processed by \code{\link{rescue.rectangle}}.
  }
  \item{\dots}{
    Optional arguments passed to \code{\link[spatstat.geom]{as.mask}}
    controlling the pixel resolution of the transformed window,
    if \code{X} is a binary pixel mask.
  }
}
\value{
  Another window (of class \code{"owin"}) representing the
  result of applying the affine transformation.
}
\details{
  The window is subjected first to the linear transformation represented by
  \code{mat} (multiplying on the left by \code{mat}),
  and then the result is translated by the vector \code{vec}.
  
  The argument \code{mat} must be a nonsingular \eqn{2 \times 2}{2 * 2}
  matrix.

  This is a method for the generic function \code{\link{affine}}.
}
\seealso{
  \code{\link{affine}},
  \code{\link{affine.ppp}},
  \code{\link{affine.psp}},
  \code{\link{affine.im}},
  \code{\link{rotate}},
  \code{\link{shift}}
}
\examples{
  # shear transformation
  shear <- matrix(c(1,0,0.6,1),ncol=2)
  X <- affine(owin(), shear)
  if(interactive()) plot(X)
  affine(letterR, shear, c(0, 0.5))
  affine(as.mask(letterR), shear, c(0, 0.5))
}
\author{\adrian
  
  
  and \rolf
  
}
\keyword{spatial}
\keyword{math}
\concept{Geometrical transformations}
