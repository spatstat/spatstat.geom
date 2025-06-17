\name{identify.tess}
\alias{identify.tess}
\title{Interactively Identify Tiles of a Tessellation}
\description{
  If a tessellation is plotted in the graphics window,
  then each time the left mouse button is pressed,
  this function will find the tile which contains
  the mouse position, and print its mark value (or its serial number
  if there is no mark).
}
\usage{
  \method{identify}{tess}(x, \dots,
                  labels=tilenames(x),
                  n=nobjects(x), plot=TRUE, paint=plot, paint.args=list())
}
\arguments{
  \item{x}{
    A tessellation
    (object of class \code{"tess"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link{plot.owin}}
    and \code{\link[graphics]{text.default}} controlling the graphical
    display.
  }
  \item{labels}{
    Labels associated with the tiles, to be plotted when the
    tiles are identified. A character vector or numeric vector
    of length equal to the number of tiles in \code{x}.
  }
  \item{n}{
    Maximum number of tiles to be identified.
  }
  \item{plot}{
    Logical. Whether to plot the labels when a tile is identified.
  }
  \item{paint}{
    Logical. Whether to redraw each identified tile, using a
    different colour.
  }
  \item{paint.args}{
    Optional list of arguments passed to \code{\link{plot.owin}} determining the
    colour and style in which each identified tile will be
    redrawn, if \code{paint=TRUE}.
  }

}
\value{
  A \code{data.frame} with columns \code{id} and \code{name}
  containing the serial numbers and names of the tiles of
  \code{x} that were identified, in the order that they were identified;
  If \code{x} is marked, subsequent columns contain the marks for these tiles.
}
\details{
  This is a method for the generic function \code{\link[graphics]{identify}}
  for tessellation objects.

  The tessellation \code{x} should first be plotted
  using \code{\link{plot.tess}}.
  Then \code{identify(x)}
  reads the position of the graphics pointer each time the
  left mouse button is pressed.  It then determines which tile of
  \code{x} contains the mouse position.
  The index of this tile (and its mark if any) 
  will be returned as part of the value of the call.

  Each time a tile is identified,
  the tile will be plotted in light blue fill colour, and 
  text will be displayed inside the tile
  showing the name of the tile.
}
\seealso{
  \code{\link[graphics]{identify}},
  \code{\link{plot.tess}}
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{iplot}
