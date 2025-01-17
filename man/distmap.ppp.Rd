\name{distmap.ppp}  
\alias{distmap.ppp}
\title{
  Distance Map of Point Pattern
}
\description{
  Computes the distance from each pixel to the nearest
  point in the given point pattern.
}
\usage{
  \method{distmap}{ppp}(X, \dots, clip=FALSE, metric=NULL)
}
\arguments{
  \item{X}{A point pattern (object of class \code{"ppp"}).
  }
  \item{\dots}{Arguments passed to \code{\link[spatstat.geom]{as.mask}}
    to control pixel resolution.
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
  The ``distance map'' of a point pattern \eqn{X} is the function
  \eqn{f} whose value \code{f(u)} is defined for any two-dimensional
  location \eqn{u} as the shortest distance from \eqn{u} to \eqn{X}.

  This function computes the distance map of the point pattern \code{X}
  and returns the distance map as a pixel image. The greyscale value
  at a pixel \eqn{u} equals the distance from \eqn{u} to the nearest
  point of the pattern \code{X}.

  If \code{clip=FALSE} (the default), the resulting pixel values are defined at every pixel in
  the rectangle \code{Frame(X)}. If \code{clip=TRUE}, the pixel values are defined
  only inside \code{Window(X)}, and are \code{NA} outside this window.
  Computation is faster when \code{clip=FALSE}.
  
  Additionally, the return value 
  has two attributes, \code{"index"} and \code{"bdry"}, which are
  also pixel images. The grey values in \code{"bdry"} give the
  distance from each pixel to the boundary of the window containing \code{X}.
  The grey values in \code{"index"} are integers identifying which
  point of \code{X} is closest. 
  
  This is a method for the generic function \code{\link{distmap}}.
  
  Note that this function gives the distance from the
  \emph{centre of each pixel} to the nearest data point.
  To compute the exact distance from a given spatial location
  to the nearest data point in \code{X}, use \code{\link{distfun}} or 
  \code{\link{nncross}}.
}
\section{Distance values}{
  The pixel values in the image \code{distmap(X)} are distances,
  expressed as multiples of the unit of length of the spatial
  coordinates in \code{X}.
  The unit of length is given by \code{\link{unitname}(X)}.

  Note that, if the unit of length in \code{X} is a composite
  expression such as \sQuote{2 microns},
  then the values in \code{distmap(X)} are expressed as multiples
  of 2 microns, rather than being expressed in microns.
}
\seealso{
  Generic function \code{\link{distmap}} and other methods
  \code{\link{distmap.psp}},
  \code{\link{distmap.owin}}.
  
  Generic function \code{\link{distfun}}.

  Nearest neighbour distance \code{\link{nncross}}.

  \code{\link{unitname}} and \code{\link{rescale}} to control the unit
  of length.
}
\examples{
  U <- distmap(cells)
  if(interactive()) {
    plot(U)
    plot(attr(U, "bdry"))
    plot(attr(U, "index"))
  }
}
\author{
  \adrian
  and \rolf
}
\keyword{spatial}
\keyword{math}
