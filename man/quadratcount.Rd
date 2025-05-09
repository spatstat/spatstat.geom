\name{quadratcount}
\alias{quadratcount}
\alias{quadratcount.ppp}
\alias{quadratcount.splitppp}
\title{Quadrat counting for a point pattern}
\description{
  Divides window into quadrats and 
  counts the numbers of points in each quadrat.
}
\usage{
  quadratcount(X, \dots)

  \method{quadratcount}{ppp}(X, nx=5, ny=nx, \dots,
               xbreaks=NULL, ybreaks=NULL, left.open=TRUE,
               tess=NULL)

  \method{quadratcount}{splitppp}(X, \dots)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"})
    or a split point pattern (object of class \code{"splitppp"}).
  }
  \item{nx,ny}{
    Numbers of rectangular quadrats in the \eqn{x} and \eqn{y} directions.
    Incompatible with \code{xbreaks} and \code{ybreaks}.
  }
  \item{\dots}{Additional arguments passed to \code{quadratcount.ppp}.}
  \item{xbreaks}{
    Numeric vector giving the \eqn{x} coordinates of the
    boundaries of the rectangular quadrats. Incompatible with \code{nx}.
  }
  \item{ybreaks}{
    Numeric vector giving the \eqn{y} coordinates of the
    boundaries of the rectangular quadrats. Incompatible with \code{ny}.
  }
  \item{tess}{
    Tessellation (object of class \code{"tess"} or something acceptable
    to \code{\link{as.tess}}) determining the quadrats. Incompatible
    with \code{nx,ny,xbreaks,ybreaks}.
  }
  \item{left.open}{
    Logical value specifying whether rectangular quadrats are left-open and
    right-closed (\code{left.open=TRUE}, the default) or
    left-closed and right-open (\code{left.open=FALSE}).
  }
}
\value{
  The value of \code{quadratcount.ppp} is a 
  contingency table containing the number of points in each
  quadrat. The table is also an object of the
  special class \code{"quadratcount"}
  and there is a plot method for this class.

  The value of \code{quadratcount.splitppp} is a list of such
  contingency tables, each containing the quadrat counts for one of the
  component point patterns in \code{X}.
  This list also has the class \code{"solist"} which has
  print and plot methods.
}
\details{
  Quadrat counting is an elementary technique for analysing spatial
  point patterns. See Diggle (2003).

  \bold{If \code{X} is a point pattern}, then 
  by default, the window containing the point pattern \code{X} is divided into
  an \code{nx * ny} grid of rectangular tiles or `quadrats'.
  (If the window is not a rectangle, then these tiles are intersected
  with the window.)
  The number of points of \code{X} falling in each quadrat is
  counted. These numbers are returned as a contingency table.

  If \code{xbreaks} is given, it should be a numeric vector
  giving the \eqn{x} coordinates of the quadrat boundaries.
  If it is not given, it defaults to a
  sequence of \code{nx+1} values equally spaced
  over the range of \eqn{x} coordinates in the window \code{Window(X)}.

  Similarly if \code{ybreaks} is given, it should be a numeric
  vector giving the \eqn{y} coordinates of the quadrat boundaries.
  It defaults to a vector of \code{ny+1} values
  equally spaced over the range of \eqn{y} coordinates in the window.
  The lengths of \code{xbreaks} and \code{ybreaks} may be different.

  Alternatively, quadrats of any shape may be used.
  The argument \code{tess} can be a tessellation (object of class
  \code{"tess"}) whose tiles will serve as the quadrats.
  
  The algorithm counts the number of points of \code{X}
  falling in each quadrat, and returns these counts as a
  contingency table. 

  The return value is a \code{table} which can be printed neatly.
  The return value is also a member of the special class
  \code{"quadratcount"}. Plotting the object will display the
  quadrats, annotated by their counts. See the examples.

  To perform a chi-squared test based on the quadrat counts,
  use \code{\link[spatstat.explore]{quadrat.test}}.
  
  To calculate an estimate of intensity based on the quadrat counts,
  use \code{\link{intensity.quadratcount}}.

  To extract the quadrats used in a \code{quadratcount} object,
  use \code{\link{as.tess}}.

  \bold{If \code{X} is a split point pattern} (object of class
  \code{"splitppp"} then quadrat counting will be performed on
  each of the components point patterns, and the resulting
  contingency tables will be returned in a list. This list can be
  printed or plotted.

  Marks attached to the points are ignored by \code{quadratcount.ppp}.
  To obtain a separate contingency table for each type of point
  in a multitype point pattern,
  first separate the different points using \code{\link{split.ppp}},
  then apply \code{quadratcount.splitppp}. See the Examples.
}
\section{Treament of data points on the boundary}{
  If the quadrats are rectangular, they are assumed to be 
  left-open and right-closed, by default (\code{left.open=TRUE}).
  Alternatively if \code{left.open=FALSE} then rectangular quadrats are
  left-closed and right-open.

  If the quadrats are not rectangular, the treatment of points which lie
  on the boundary of two quadrats is undefined, and may depend on the
  hardware.
}
\section{Hypothesis testing}{
  To perform a chi-squared test based on the quadrat counts,
  use \code{\link[spatstat.explore]{quadrat.test}}.
}
\section{Warning}{
  If \code{Q} is a \code{quadratcount} object,
  the ordering of entries in the table \code{Q} 
  \bold{may be different from} the ordering of quadrats (tiles
  in the tessellation \code{as.tess(Q)}).

  To obtain the entries of the table in the same order
  as the quadrats, use
  \code{counts <- as.numeric(t(Q))} or \code{counts <- marks(as.tess(Q))}.
}
\seealso{
  \code{\link{plot.quadratcount}},
  \code{\link{intensity.quadratcount}},
  \code{\link{quadrats}},
  \code{\link[spatstat.explore]{quadrat.test}},
  \code{\link{tess}},
  \code{\link{hextess}},
  \code{\link[spatstat.random]{quadratresample}},
  \code{\link[spatstat.explore]{miplot}}
}
\references{
  Diggle, P.J. \emph{Statistical analysis of spatial point patterns}.
  Academic Press, 2003.

  Stoyan, D. and Stoyan, H. (1994)
  \emph{Fractals, random shapes and point fields:
  methods of geometrical statistics.}
  John Wiley and Sons.
} 
\examples{
 X <- runifrect(50)
 quadratcount(X)
 quadratcount(X, 4, 5)
 quadratcount(X, xbreaks=c(0, 0.3, 1), ybreaks=c(0, 0.4, 0.8, 1))
 qX <-  quadratcount(X, 4, 5)

 # plotting:
 plot(X, pch="+")
 plot(qX, add=TRUE, col="red", cex=1.5, lty=2)

 # irregular window
 plot(humberside)
 qH <- quadratcount(humberside, 2, 3)
 plot(qH, add=TRUE, col="blue", cex=1.5, lwd=2)

 # multitype - split
 plot(quadratcount(split(humberside), 2, 3))
 
 # quadrats determined by tessellation:
 B <- dirichlet(runifrect(6))
 qX <- quadratcount(X, tess=B)
 plot(X, pch="+")
 plot(qX, add=TRUE, col="red", cex=1.5, lty=2)
}
\author{\adrian
  and \rolf
}
\keyword{spatial}
\keyword{math}
 
 
