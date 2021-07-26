\name{convexpairdist}
\alias{convexpairdist}
\alias{convexcrossdist}
\alias{convexnndist}
\alias{convexnnwhich}
\alias{convexnncross}
\title{
  Distances Using Arbitrary Convex Metric
}
\description{
  Computes the pairwise distances or nearest-neighbour distances
  between points in a point pattern, using the distance
  metric associated with a convex polygon.
}
\usage{
convexpairdist(X, K)
convexnndist(X, K)
convexnnwhich(X, K)
convexcrossdist(X, Y, K)
convexnncross(X, Y, K, what=c("dist", "which"))
}
\arguments{
  \item{X,Y}{
    Point patterns (objects of class \code{"ppp"}).
  }
  \item{K}{
    Convex set defining the metric.
    A polygon that is symmetric about the origin. See Details.
  }
  \item{what}{
    Character string or character vector (partially matched,
    several entries accepted) specifying which values to compute.
  }
}
\details{
  These functions compute distances between points in two-dimensional
  space using the distance metric whose unit ball is \code{K}.
  \itemize{
    \item \code{convexpairdist} is the counterpart of
    \code{\link{pairdist}}. It computes the distances between all pairs
    of points in the pattern \code{X}, and returns a square matrix of
    distance values.
    \item \code{convexcrossdist} is the counterpart of
    \code{\link{crossdist}}.
    It computes the distances from each point in the pattern \code{X}
    to each point in the pattern \code{Y}, and returns a matrix of
    distance values.
    \item \code{convexnndist} is the counterpart of
    \code{\link{nndist}}.
    It computes the distance from each point of \code{X} to the nearest
    other point of \code{X}, and returns a vector of distance values.
    \item \code{convexnnwhich} is the counterpart of \code{\link{nnwhich}}.
    For each point of \code{X}, it finds the nearest other point of
    \code{X}, and returns a vector of integer indices.
    \item \code{convexnncross} is the counterpart of
    \code{\link{nncross}}. For each point of \code{X}, it finds the
    distance to the nearest point of \code{Y}, and/or the index of the
    nearest point of \code{Y}. It returns a vector or a two-column data
    frame.
  }
  \code{K} will be converted to a polygon.
  It must be convex, and symmetric about the origin.
}
\value{
  A matrix, vector or data frame of numeric values.
}
\author{
  \adrian
}
\seealso{
  \code{\link{pairdist}},
  \code{\link{crossdist}},
  \code{\link{nndist}},
  \code{\link{nnwhich}},
  \code{\link{nncross}}.

  \code{\link{convexdistmap}}
}
\examples{
  K <- owin(poly=list(x=c(2.5,2,0.5,-2.5,-2,-0.5),y=c(0,1,2,0,-1,-2)))
  plot(K)
  plot(redwood)
  ## show point number 43
  plot(redwood[43], pch=16, add=TRUE)
  ## nearest neighbour distances
  nn <- convexnndist(redwood, K)
  ## Nearest neighbour distance for point number 43 is nn[43]; verify
  plot(shift(scalardilate(K, nn[43]), redwood[43]), add=TRUE)
  kk <- convexnnwhich(redwood, K)
  ## nearest neighbour for point number 43 is point number kk[43]; verify
  plot(redwood[kk[43]], pch=3, col="red", add=TRUE)
}
\keyword{spatial}
\keyword{math}
