\name{bufftess}
\alias{bufftess}
\title{
  Buffer Distance Tessellation
}
\description{
  Constructs a spatial tessellation, composed of rings or buffers at specified
  distances away from the given spatial object.
}
\usage{
bufftess(X, breaks, W = Window(X), \dots, polygonal = TRUE)
}
\arguments{
  \item{X}{
    A spatial object in two dimensions,
    such as a point pattern (class \code{"ppp"}) or
    line segment pattern (class \code{"psp"}).
  }
  \item{breaks}{
    Either a numeric vector specifying the cut points for the
    distance values, or a single integer specifying the number of cut points.
  }
  \item{W}{
    Optional. Window (object of class \code{"owin"}) inside which the
    tessellation will be constructed.
  }
  \item{\dots}{
    Optional arguments passed to \code{\link{as.mask}} controlling the
    pixel resolution when \code{polygonal=FALSE},
    and optional arguments passed to \code{\link[base]{cut.default}}
    controlling the labelling of the distance bands.
}
  \item{polygonal}{
    Logical value specifying whether the tessellation should consist of
    polygonal tiles (\code{polygonal=TRUE}, the default) or
    should be constructed using a pixel image (\code{polygonal=FALSE}).
  }
}
\details{
  This function divides space into tiles defined
  by distance from the object \code{X}. The result is a tessellation
  (object of class \code{"tess"}) that consists of concentric rings
  around \code{X}.

  The distance values which determine the tiles are specified by
  the argument \code{breaks}.
  \itemize{
    \item
    If \code{breaks} is a vector of numerical values, then these
    values are taken to be the distances defining the tiles.
    The first tile is the region of space that lies at distances
    between \code{breaks[1]} and \code{breaks[2]} away from \code{X};
    the second tile is the region lying at distances between
    \code{breaks[2]} and \code{breaks[3]} away from \code{X}; and so on.
    The number of tiles will be \code{length(breaks)-1}.
    \item
    If \code{breaks} is a single integer, it is interpreted as specifying
    the number of intervals between breakpoints.
    There will be \code{breaks+1} equally spaced
    break points, ranging from zero to the maximum
    achievable distance. The number of tiles will equal \code{breaks}.
  }

  The tessellation can be computed using either raster calculations
  or vector calculations.
  \itemize{
    \item 
    If \code{polygonal=TRUE} (the default), the tiles are computed as
    polygonal windows using vector geometry, and the result is a
    tessellation consisting of polygonal tiles.
    This calculation could be slow and could require substantial memory,
    but produces a geometrically accurate result.
    \item
    If \code{polygonal=FALSE}, the distance map of \code{X} is computed
    as a pixel image (\code{\link{distmap}}), then the distance values
    are divided into discrete bands using \code{\link{cut.im}}. The result is a
    tessellation specified by a pixel image. This computation is faster
    but less accurate.
  }
}
\value{
  A tessellation (object of class \code{"tess"}).
  
  The result also has an attribute \code{breaks} which is the vector of distance breakpoints.
}
\author{
  \spatstatAuthors.
}
\seealso{
  Polygonal calculations are performed using \code{\link{dilation}}
  and \code{\link{setminus.owin}}. Pixel calculations are performed
  using \code{\link{distmap}} and \code{\link{cut.im}}. See
  \code{\link{as.mask}} for details of arguments that control pixel
  resolution. 
  
  For other kinds of tessellations, see 
  \code{\link{tess}},
  \code{\link{hextess}},
  \code{\link{venn.tess}},
  \code{\link{polartess}},
  \code{\link{dirichlet}}, \code{\link{delaunay}},
  \code{\link{quantess}}, \code{\link{quadrats}}
  and 
  \code{\link[spatstat.random]{rpoislinetess}}.
}
\examples{
  X <- cells[c(FALSE,FALSE,FALSE,TRUE)] 
  b <- c(0, 0.05, 0.1, 0.15, 0.2, Inf)
  plot(bufftess(X, b), do.col=TRUE, col=1:5)
}
\keyword{spatial}
\keyword{manip}
