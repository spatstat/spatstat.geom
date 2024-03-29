\name{plot.quadratcount}
\alias{plot.quadratcount}
\title{
  Plot Quadrat Counts
}
\description{
  Given a table of quadrat counts for a spatial point pattern,
  plot the quadrats which were used, and display the quadrat count
  as text in the centre of each quadrat.
}
\usage{
\method{plot}{quadratcount}(x, \dots, add = FALSE,
              entries=as.integer(t(x)), 
              dx = 0, dy = 0, show.tiles = TRUE,
              textargs = list())
}
\arguments{
  \item{x}{
    Object of class \code{"quadratcount"}
    produced by the function \code{\link{quadratcount}}.
  }
  \item{\dots}{
    Additional arguments passed to \code{\link{plot.tess}}
    to plot the quadrats.
  }
  \item{add}{
    Logical. Whether to add the graphics to an existing plot.
  }
  \item{entries}{
    Vector of numbers to be plotted in each quadrat.
    The default is to plot the quadrat counts.
  }
  \item{dx,dy}{
    Horizontal and vertical displacement of text
    relative to centroid of quadrat.
  }
  \item{show.tiles}{
    Logical value indicating whether to plot the quadrats.
  }
  \item{textargs}{
    List containing extra arguments
    passed to \code{\link[graphics]{text.default}}
    to control the annotation.
  }
}
\details{
  This is the plot method for the objects 
  of class \code{"quadratcount"} that are
  produced by the function \code{\link{quadratcount}}.
  Given a spatial point pattern, \code{\link{quadratcount}}
  divides the observation window into disjoint tiles or quadrats,
  counts the number of points in each quadrat, and stores the
  result as a contingency table which also belongs to the
  class \code{"quadratcount"}.

  First the quadrats are plotted
  (provided \code{show.tiles=TRUE}, the default).
  This display can be controlled by passing additional arguments \code{\dots}
  to \code{\link{plot.tess}}.
  
  Then the quadrat counts are printed using
  \code{\link[graphics]{text.default}}. This display can be controlled
  using the arguments \code{dx,dy} and \code{textargs}.

  If \code{entries} is given, it should be a vector of length equal to
  the number of quadrats (the number of tiles in the tessellation
  \code{as.tess(x)}) containing integer or character values
  to be displayed in each quadrat, in the same sequence
  as \code{tiles(as.tess(x))} or \code{tilenames(as.tess(x))}
  or the counts in the transposed table \code{t(x)}.
}
\value{
  Null.
}
\seealso{
  \code{\link{quadratcount}},
  \code{\link{plot.tess}},
  \code{\link[graphics]{text.default}},
  \code{\link[spatstat.explore]{plot.quadrattest}}.
}
\examples{
   plot(quadratcount(swedishpines, 5))
}
\author{
  \adrian
  and \rolf
}
\keyword{spatial}
\keyword{hplot}
