\name{clickpoly}
\alias{clickpoly}
\title{Interactively Define a Polygon}
\description{
  Allows the user to create a polygon by
  point-and-click in the display.
}
\usage{
  clickpoly(add=FALSE, nv=NULL, np=1, \dots, snap.step=NULL, snap.origin=c(0,0))
}
\arguments{
  \item{add}{
    Logical value indicating whether to create a new plot
    (\code{add=FALSE}) or draw over the existing plot (\code{add=TRUE}).
  }
  \item{nv}{
    Number of vertices of the polygon (if this is predetermined).
  }
  \item{np}{
    Number of polygons to create.
  }
  \item{\dots}{
    Arguments passed to \code{\link[graphics]{locator}} to control the
    interactive plot, and to \code{\link[graphics]{polygon}} to plot the
    polygons.
  }
  \item{snap.step}{
    Optional. Spatial coordinates will be rounded to the nearest
    multiple of \code{snap.step}. A positive number specifying the step
    length, or a vector of 2 positive numbers specifying step lengths
    for the \eqn{x} and \eqn{y} coordinates.
  }
  \item{snap.origin}{
    Optional. Numeric vector of length 2. Coordinates of the origin
    that will be used when rounding coordinates.
  }
}
\value{
  A window (object of class \code{"owin"}) representing the polygon.
}
\details{
  This function allows the user to create a polygonal window 
  by interactively clicking on the screen display.

  The user is prompted to point the mouse at any desired locations
  for the polygon vertices,
  and click the left mouse button to add each point.
  Interactive input stops after \code{nv} clicks (if \code{nv} was given)
  or when the middle mouse button is pressed.

  The return value is a window (object of class \code{"owin"})
  representing the polygon.

  This function uses the \R command \code{\link[graphics]{locator}} to
  input the mouse clicks. It only works on screen devices such as
  \sQuote{X11}, \sQuote{windows} and \sQuote{quartz}. Arguments that can be
  passed to \code{\link[graphics]{locator}} through \code{\dots} include
  \code{pch} (plotting character), \code{cex} (character expansion
  factor) and \code{col} (colour). See \code{\link[graphics]{locator}}
  and \code{\link[graphics]{par}}.

  Multiple polygons can also be drawn, by specifying
  \code{np > 1}. The polygons must be disjoint. The result is 
  a single window object consisting of all the polygons.

  If \code{snap.step} is given, a grid of lines will be drawn with
  the specified spacing, and the clicked locations will be 
  snapped (rounded) to the nearest grid point.
}
\seealso{
  \code{\link{identify.ppp}},
  \code{\link{clickbox}},
  \code{\link{clickppp}},
  \code{\link{clickdist}},
  \code{\link[graphics]{locator}}
}
\author{
  \adrian
  and \rolf.
}
\keyword{spatial}
\keyword{iplot}
