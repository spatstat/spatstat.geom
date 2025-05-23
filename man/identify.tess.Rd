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
  \method{identify}{tess}(x, \dots)
}
\arguments{
  \item{x}{
    A tessellation
    (object of class \code{"tess"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link[graphics]{identify.default}}.
  }
}
\value{
  If \code{x} is unmarked, the result is 
  a vector containing the serial numbers of the tiles of
  \code{x} that were identified.
  If \code{x} is marked, the result is a 
  data frame, the first column containing the serial numbers
  and the subsequent columns containing the marks for these tiles.
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
  text will be displayed inside the tile
  showing its serial number.
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
