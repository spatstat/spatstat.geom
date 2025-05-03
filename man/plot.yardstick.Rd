\name{plot.yardstick}
\alias{plot.yardstick}
\title{
  Plot a Yardstick or Scale Bar
}
\description{
  Plots an object of class \code{"yardstick"}.
}
\usage{
 \method{plot}{yardstick}(x, \dots,
        style=c("arrows", "zebra"),
        angle = 20, frac = 1/8, 
        split = FALSE, shrink = 1/4,
        pos = NULL,
        txt.args=list(),
        txt.shift=c(0,0),
        zebra.step=NULL, zebra.width=NULL,
        zebra.col="black", zebra.scale=1,
        zebra.args=list(), zebra.shift=c(0,0),
        do.plot = TRUE, do.txt=TRUE)
}
\arguments{
  \item{x}{
    Object of class \code{"yardstick"} to be plotted.
    This object is created by the command \code{\link{yardstick}}.
  }
  \item{\dots}{
    Additional graphics arguments passed to
    \code{\link[graphics]{segments}} to control the appearance of the line.
  }
  \item{style}{
    Character string (partially matched) specifying the
    style of plot. See Details.
  }
  \item{angle}{
    Angle between the arrows and the line segment, in degrees.
    Applies when \code{style="arrows"}.
  }
  \item{frac}{
    Length of arrow as a fraction of total length of the line segment.
    Applies when \code{style="arrows"}.
  }
  \item{split}{
    Logical. If \code{TRUE}, then the line will be broken in the
    middle, and the text will be placed in this gap. If \code{FALSE},
    the line will be unbroken, and the text will be placed beside
    the line.
    Applies when \code{style="arrows"}.
  }
  \item{shrink}{
    Fraction of total length to be removed from the middle of the
    line segment, if \code{split=TRUE}.
    Applies when \code{style="arrows"}.
  }
  \item{pos}{
    Integer (passed to \code{\link[graphics]{text}}) determining the
    position of the annotation text relative to the line segment,
    if \code{split=FALSE}. Values of 1, 2, 3 and 4
    indicate positions below, to the left of,
    above and to the right of the line, respectively.
  }
  \item{txt.args}{
    Optional list of additional arguments passed to
    \code{\link[graphics]{text}} controlling the appearance of the text.
    Examples include \code{adj}, \code{srt}, \code{col}, \code{cex},
    \code{font}.
  }
  \item{txt.shift}{
    Optional numeric vector of length 2 specifying displacement
    of the text position relative to the centre of the yardstick.
  }
  \item{zebra.step}{
    Length of each bar in the zebra pattern.
    Applies when \code{style="zebra"}.
  }
  \item{zebra.width}{
    Width of each bar in the zebra pattern.
    Applies when \code{style="zebra"}.
  }
  \item{zebra.col}{
    Colour of each bar in the zebra pattern.
    Applies when \code{style="zebra"}.
  }
  \item{zebra.scale}{
    Scale value for numerical labels in the zebra pattern.
    Physical lengths will be divided by \code{zebra.scale}
    to determine the corresponding numerical labels.
  }
  \item{zebra.args}{
    Optional list of additional arguments passed to
    \code{\link[graphics]{text}} controlling the appearance of the
    text annotation for each bar in the zebra pattern.
    Examples include \code{adj}, \code{srt}, \code{col}, \code{cex},
    \code{font}.
  }
  \item{zebra.shift}{
    Optional numeric vector of length 2 specifying displacement
    of the text annotation for each bar in the zebra pattern.
  }
  \item{do.plot}{
    Logical values specifying whether
    to actually perform the plot.
  }
  \item{do.txt}{
    Logical value specifying whether to draw text annotation.
  }
}
\details{
  A yardstick or scale bar is a line segment, drawn on any spatial
  graphics display, indicating the scale of the plot.

  \itemize{
    \item If \code{style="arrows"}, the line segment is drawn as a pair
    of arrows pointing from the middle of the line to the ends of the
    line. This style is often used in architectural drawings.
    If \code{angle=0}, the arrow heads are replaced by parallel bars
    marking the two ends of the line.
    \item If \code{style="zebra"}, the line segment is divided into
    block of length \code{zebra.step} and width \code{zebra.width}
    units. Blocks are drawn alternately as filled rectangles and
    outlined rectangles, so that the result resembles a zebra crossing.
    This style is often used in maps and charts. There are sensible
    defaults for \code{zebra.step} and \code{zebra.width}.
  }
  
  The argument \code{x} should be an object of class \code{"yardstick"}
  created by the command \code{\link{yardstick}}.

}
\value{
  A window (class \code{"owin"}) enclosing the plotted graphics.
}
\examples{
  plot(owin(), main="Yardsticks")
  ys <- yardstick(as.psp(list(xmid=0.5, ymid=0.1, length=0.4, angle=0),
                         window=owin(c(0.2, 0.8), c(0, 0.2))),
                  txt="1 km")
  plot(ys)
  ys <- shift(ys, c(0, 0.3))
  plot(ys, angle=90, frac=0.08)
  ys <- shift(ys, c(0, 0.3))
  plot(ys, split=TRUE)

  yt <- shift(ys, c(0, 0.2))
  plot(yt, style="z", pos=3, zebra.step=0.1, txt.args=list(offset=0.1))
}
\author{\spatstatAuthors.}
\seealso{
  \code{\link{yardstick}}
}
\keyword{spatial}
\keyword{hplot}

