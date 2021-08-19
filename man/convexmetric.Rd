\name{convexmetric}
\alias{convexmetric}
\title{
  Distance Metric Defined by Convex Set
}
\description{
  Create the distance metric associated with a given convex polygon.
}
\usage{
convexmetric(K)
}
\arguments{
  \item{K}{
    Convex set defining the metric.
    A polygon that is symmetric about the origin. See Details.
  }
}
\details{
  This function creates the distance metric
  associated with the convex set \code{K}
  so that the unit ball of the metric is equal to \code{K}.
  It returns an object of class \code{"metric"} representing the metric
  (see \code{\link{metric.object}}).
  
  The argument \code{K} must be a window (class \code{"owin"}).
  It will be converted to a polygon.
  It must be convex, and symmetric about the origin.
}
\value{
  An object of class \code{"metric"}.
}
\author{
  \adrian
}
\seealso{
  \code{\link{metric.object}}
}
\examples{
  K <- owin(poly=list(x=c(2.5,2,0.5,-2.5,-2,-0.5),y=c(0,1,2,0,-1,-2)))
  plot(K)
  points(0,0)
  m <- convexmetric(K)
  ## show redwood data and identify point number 43
  plot(redwood, main="")
  plot(redwood[43], pch=16, add=TRUE)
  ## nearest neighbour distances and identifiers
  nd <- m[["nndist"]](redwood)
  nw <- m[["nnwhich"]](redwood)
  ## Nearest neighbour distance for point number 43 is nd[43]; verify
  B43 <- disc(radius=nd[43], centre=redwood[43], metric=m)
  plot(B43, add=TRUE)
  ## nearest neighbour for point number 43 is point number nw[43]; verify
  plot(redwood[nw[43]], pch=3, col="red", add=TRUE)
}
\keyword{spatial}
\keyword{math}
