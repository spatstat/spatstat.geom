\name{rounding.ppp} 
\alias{rounding.ppp}
\alias{rounding.pp3}
\alias{rounding.ppx}
\title{
  Detect Numerical Rounding 
}
\description{
  Given a numeric vector, or an object containing numeric spatial
  coordinates, determine whether the values have been rounded
  to a certain number of decimal places.
}
\usage{
\method{rounding}{ppp}(x)

\method{rounding}{pp3}(x)

\method{rounding}{ppx}(x)
}
\arguments{
  \item{x}{
    A point pattern (object of class \code{ppp}, \code{pp3} or \code{ppx}).
  }
}
\details{
  The functions documented here are methods for the generic
  \code{\link[spatstat.univar]{rounding}}.
  They determine whether the coordinates of a spatial object
  have been rounded to a certain number of decimal places.
  \itemize{
    \item 
    If the coordinates of the points in \code{x} are not all integers, then
    \code{rounding(x)} returns the smallest number of digits \code{d}
    after the decimal point
    such that \code{\link[base]{round}(coords(x), digits=d)} is identical to
    \code{coords(x)}.
    For example if \code{rounding(x) = 2} then the coordinates of the
    points in \code{x} appear to have been rounded to 2 decimal places,
    and are multiples of 0.01.
    \item 
    If all the coordinates of the points in \code{x} are integers, then
    \code{rounding(x)} returns \code{-d}, where
    \code{d} is the smallest number of digits \emph{before} the decimal point
    such that \code{\link[base]{round}(coords(x), digits=-d)} is identical to
    \code{coords(x)}.
    For example if \code{rounding(x) = -3} then the coordinates of all
    points in \code{x} are multiples of 1000.
    If \code{rounding(x) = 0} then the entries of \code{x} are integers
    but not multiples of 10.
    \item
    If all coordinates of points in \code{x} are equal to 0,
    a value of 0 is returned.
  }
}
\value{
  An integer.
}
\author{
  \adrian
  and \rolf
}
\seealso{
  \code{\link{round.ppp}},
  \code{\link[spatstat.univar]{rounding}}
}
\examples{
   rounding(cells)
}
\keyword{spatial}
\keyword{math}
