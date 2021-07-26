\name{convexdistmap}
\alias{convexdistmap}
\title{
  Distance Map Using Arbitrary Convex Metric
}
\description{
  Computes the distance map of a spatial region
  based on the metric associated with a convex polygon.
}
\usage{
convexdistmap(W, K, npasses=5, verbose=FALSE)
}
\arguments{
  \item{W}{
    A window (object of class \code{"owin"}).
  }
  \item{K}{
    Convex set defining the metric.
    A polygon that is symmetric about the origin. See Details.
  }
  \item{npasses}{
    Experimental. 
  }
  \item{verbose}{
    Logical value indicating whether to print trace information.
  }
}
\details{
  This function computes the distance map of the
  spatial region \code{W} using the distance metric
  whose unit ball is \code{K}.

  \code{K} will be converted to a polygon.
  It must be convex, and symmetric about the origin.
}
\value{
  A pixel image (object of class \code{"im"}).
}
\author{
  \adrian
}
\seealso{
  \code{\link{distmap}}, \code{\link{rectdistmap}}
}
\examples{
  V <- letterR
  Frame(V) <- grow.rectangle(Frame(V), 0.5)
  K <- owin(poly=list(x=c(2.5,2,0.5,-2.5,-2,-0.5),y=c(0,1,2,0,-1,-2)))
  plot(K)
  plot(convexdistmap(V, K))
}
\keyword{spatial}
\keyword{math}
