\name{persp.ppp}
\alias{persp.ppp}
\title{
  Perspective Plot of Marked Point Pattern
}
\description{
  For a spatial point pattern with numeric marks,
  generate a perspective plot in which each data point is
  shown as a vertical spike, with height proportional to the mark value.
}
\usage{
  \method{persp}{ppp}(x, \dots, main, grid = TRUE, ngrid = 10,
            col.grid = "grey", col.base = "white",
            win.args=list(), spike.args = list(), neg.args = list(),
            which.marks = 1, zlab = NULL, zlim = NULL, zadjust = 1,
            legend=TRUE, legendpos="bottomleft",
            leg.args=list(lwd=4), leg.col=c("black", "orange"))
}
\arguments{
  \item{x}{
    A spatial point pattern (object of class \code{"ppp"})
    with numeric marks, or a data frame of marks.
  }
  \item{\dots}{
    Additional graphical arguments passed to
    \code{\link{persp.default}} to determine the perspective view
    (for example the rotation angle \code{theta} and the elevation angle
    \code{phi})
    or passed to \code{\link[graphics]{segments}} to control the drawing
    (for example \code{lwd} for line width).
  }
  \item{main}{
    Optional main title for the plot.
  }
  \item{grid}{
    Logical value specifying whether to draw a grid of reference lines
    on the horizontal plane.
  }
  \item{ngrid}{
    Number of grid lines to draw in each direction, if \code{grid=TRUE}.
    An integer, or a pair of integers specifying the number of grid
    lines along the horizontal and vertical axes respectively.
  }
  \item{col.grid}{
    Colour of grid lines, if \code{grid=TRUE}.
  }
  \item{col.base}{
    Colour with which to fill the horizontal plane.
  }
  \item{win.args}{
    List of arguments passed to \code{\link{plot.owin}}
    to control the drawing of the window of \code{x}.
    Applicable only when the window is not a rectangle.
  }
  \item{spike.args}{
    List of arguments passed to \code{\link[graphics]{segments}}
    to control the drawing of the spikes.
  }
  \item{neg.args}{
    List of arguments passed to \code{\link[graphics]{segments}}
    applicable only to those spikes which have negative height
    (corresponding to a mark value which is negative).
  }
  \item{which.marks}{
    Integer, or character name, identifying the column of marks which
    should be used, when \code{marks(x)} is a data frame.
  }
  \item{zlab}{
    Optional. Label for the vertical axis. 
    Character string or expression.
  }
  \item{zlim}{
    Optional. Range of values on the vertical axis.
    A numeric vector of length 2.
  }
  \item{zadjust}{
    Scale adjustment factor controlling the height of spikes.
  }
  \item{legend}{
    Logical value specifying whether to draw a reference scale bar for the
    vertical axis.
  }
  \item{legendpos}{
    Position of the reference scale bar. Either a character string
    matching one of the options
    \code{"bottomleft"}, \code{"bottomright"},
    \code{"topleft"}, \code{"topright"},
    \code{"bottom"}, \code{"left"}, \code{"top"} or \code{"right"},
    or a numeric vector of length 2 specifing the coordinate position
    of the base of the reference scale bar.
  }
  \item{leg.args}{
    Additional arguments passed to \code{\link[graphics]{segments}}
    to control the drawing of the reference scale bar.
  }
  \item{leg.col}{
    A vector (usually of length 2) of colour values
    for successive intervals in the reference scale.
    The default is a reference scale consisting of
    black and orange stripes.
  }
}
\details{
  The function \code{\link[graphics]{persp}} is generic. This is the
  method for spatial point patterns (objects of class \code{"ppp"}).
  The argument \code{x} must be a point pattern with numeric marks,
  or with a data frame of marks.
  
  A perspective view will be plotted. The eye position is determined by
  the arguments \code{theta} and \code{phi}
  passed to \code{\link[graphics]{persp.default}}.

  First the horizontal plane is drawn in perspective view,
  using a faint grid of lines to help suggest the perspective.
  Next the observation window of \code{x} is placed on the horizontal
  plane and its edges are drawn in perspective view.
  Finally for each data point in \code{x}, a vertical spike is erected
  at the spatial location of the data point, with height equal to the
  mark value of the point.

  If any mark values are negative, the corresponding spikes will
  penetrate below the horizontal plane. They can be drawn in a different
  colour by specifying \code{neg.args} as shown in the examples.

  Like all spatial plots in the \pkg{spatstat} family,
  \code{persp.ppp} does not independently rescale
  the \eqn{x} and \eqn{y} coordinates. A long narrow window will be
  represented as a long narrow window in the perspective view.
  To override this and allow the coordinates to be independently
  rescaled, use the argument \code{scale=TRUE} which will be passed
  to \code{\link[graphics]{persp.default}}.
}
\value{
  (Invisibly) the perspective transformation matrix.
}
\author{
Adrian Baddeley.
}
\examples{
   persp(longleaf, theta=-30, phi=35, spike.args=list(lwd=3), zadjust=1.5)

   # negative mark values
   X <- longleaf
   marks(X) <- marks(X) - 20
   persp(X, theta=80, phi=35, neg.args=list(col="red"),
         spike.args=list(lwd=3), zadjust=1.2)

   # irregular window
   Australia <- Window(austates)
   Y <- runifrect(70, Frame(Australia))[Australia]
   marks(Y) <- runif(npoints(Y))
   persp(Y, theta=30, phi=20, col.base="lightblue",
         win.args=list(col="pink", border=NA),
         spike.args=list(lwd=2), zadjust=1.5)

}
\keyword{hplot}
\keyword{spatial}
