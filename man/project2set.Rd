\name{project2set}
\alias{project2set}
\title{
  Find Nearest Point in a Region
}
\description{
  For each data point in a point pattern \code{X},
  find the nearest location in a given spatial region \code{W}.
}
\usage{
  project2set(X, W, \dots, polygonal=is.polygonal(W))
}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"}).
  }
  \item{W}{
    Window (object of class \code{"owin"})
    or something acceptable to \code{\link{as.owin}}.
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat.geom]{as.mask}} controlling the
    pixel resolution, if \code{polygonal=FALSE}.
  }
  \item{polygonal}{
    Logical value specifying whether \code{W} should be
    treated as a polygonal window (\code{polygonal=TRUE}) or
    a binary mask window (\code{polygonal=FALSE}), and converted
    to the desired kind of window if necessary.
  }
}
\details{
  This function alters the point pattern \code{X} so that every point
  lies inside the window \code{W}.

  \itemize{
    \item
    If \code{polygonal=FALSE}
    the window \code{W} is first discretised as a binary mask
    using \code{\link[spatstat.geom]{as.mask}} with pixel resolution determined
    by the arguments \code{\dots}.
    Then for each data point \code{X[i]} in the point pattern \code{X},
    the algorithm finds the nearest pixel centre inside \code{W},
    and moves \code{X[i]} to this location.
    \item
    If \code{polygonal=TRUE}, 
    the window \code{W} is first converted to
    a polygonal window using \code{\link[spatstat.geom]{as.polygonal}}.
    Then for each data point \code{X[i]} in the point pattern \code{X},
    the algorithm finds the nearest location along the boundary of \code{W},
    moves \code{X[i]} to that location, and then
    perturbs \code{X[i]} slightly so that it lies
    a short distance inside the window \code{W}.
  }

  The result is a point pattern \code{Y} containing these nearest points,
  that is, \code{Y[i]} is the nearest point in \code{W} to the
  point \code{X[i]}.
}
\value{
  A point pattern (object of class \code{"ppp"}) with the same
  number of points as \code{X} in the window \code{W}.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{project2segment}},
  \code{\link{nncross}}
}
\examples{
## discrete case
  He <- heather$fine[owin(c(2.8, 7.4), c(4.0, 7.8))]
  plot(He, main="project2set")
  W <- erosion(complement.owin(He), 0.2)
if(require(spatstat.random)) {
  X <- runifpoint(4, W)
} else {
  X <- ppp(c(6.1, 4.3, 5.7, 4.7), c(5.0, 6.6, 7.5, 4.9), window=W)
}
  points(X, col="red")
  Y <- project2set(X, He)
  points(Y, col="green")
  arrows(X$x, X$y, Y$x, Y$y, angle=15, length=0.2)

## polygonal case
  A <- runifrect(100, Frame(letterR))
  A <- A[!inside.owin(A, , letterR)]
  plot(A, cols="red")
  B <- project2set(A, letterR)
  points(B, col="green")
  arrows(A$x, A$y, B$x, B$y, angle=15, length=0.2)
}
\keyword{spatial}
\keyword{math}

