\name{plot.colourmap}
\alias{plot.colourmap}
\title{Plot a Colour Map}
\description{
  Displays a colour map as a colour ribbon
}
\usage{
\method{plot}{colourmap}(x, ...,
       main, xlim = NULL, ylim = NULL, vertical = FALSE,
       axis = TRUE, side = if(vertical) "right" else "bottom",
       labelmap=NULL, gap=0.25, add=FALSE, increasing=NULL, nticks=5,
       at=NULL, box=NULL)
}
\arguments{
  \item{x}{Colour map to be plotted. An object of class \code{"colourmap"}.}
  \item{\dots}{
    Graphical arguments passed to \code{\link{image.default}}
    or \code{\link{axis}}.
  }
  \item{main}{Main title for plot. A character string.}
  \item{xlim}{
    Optional range of \code{x} values for the location of the
    colour ribbon.
  }
  \item{ylim}{
    Optional range of \code{y} values for the location of the
    colour ribbon.
  }
  \item{vertical}{
    Logical flag determining whether the colour ribbon
    is plotted as a horizontal strip (\code{FALSE}) or a vertical strip
    (\code{TRUE}).
  }
  \item{axis}{
    Logical value determining whether an axis should be plotted
    showing the numerical values that are mapped to the colours.
  }
  \item{side}{
    One of the character strings
    \code{"bottom"}, \code{"left"}, \code{"top"}  or \code{"right"},
    or one of the integers from 1 to 4,
    specifying the position of the axis tick marks and labels,
    if \code{axis=TRUE}.
  }
  \item{labelmap}{
    Function. If this is present, then the labels on the plot,
    which indicate the input values corresponding to particular colours,
    will be transformed by \code{labelmap} before being displayed
    on the plot. Typically used to simplify or shorten the
    labels on the plot.
  }
  \item{gap}{
    Distance between separate blocks of colour, as a fraction of
    the width of one block, if the colourmap is discrete.
  }
  \item{add}{
    Logical value indicating whether to add the colourmap to the
    existing plot (\code{add=TRUE}), or to start a new plot
    (\code{add=FALSE}, the default).
  }
  \item{increasing}{
    Logical value indicating whether to display the colour map
    in increasing order. See Details.
  }
  \item{nticks}{
    Optional. Integer specifying the approximate number of tick marks
    (representing different values of the numerical input) that should be
    drawn next to the colour map.
    Applies only when the colour map inputs are numeric values.
  }
  \item{at}{
    Optional. Numeric values at which tick marks should be drawn.
    Applies only when the colour map inputs are numeric values.
  }
  \item{box}{
    Optional. Logical value specifying whether to draw a black box
    around the colour ribbon. Default is \code{box=FALSE} when plotting
    separate blocks of colour (i.e. when the colourmap is discrete
    and \code{gap > 0}) and \code{box=TRUE} otherwise.
  }
}
\details{
  This is the plot method for the class \code{"colourmap"}.
  An object of this class
  (created by the function \code{\link{colourmap}})
  represents a colour map or
  colour lookup table associating colours with each data value.

  The command \code{plot.colourmap} displays the colour map as a colour
  ribbon or as a colour legend (a sequence of blocks of colour).
  This plot can be useful on its own to inspect the colour map.

  If the domain of the colourmap is an interval of real numbers,
  the colourmap is displayed as a continuous ribbon of colour.
  If the domain of the colourmap is a finite set of inputs,
  the colours are displayed as separate blocks of colour.
  The separation between blocks is equal to \code{gap} times
  the width of one block. 

  To annotate an existing plot with an explanatory colour ribbon
  or colour legend,
  specify \code{add=TRUE} and use the arguments \code{xlim}
  and/or \code{ylim} to control the physical position of the ribbon
  on the plot.

  Labels explaining the colour map are
  drawn by \code{\link[graphics]{axis}} and can be modified by
  specifying arguments that will be passed to this function.

  The argument \code{increasing} indicates whether the
  colourmap should be displayed so that the
  input values are increasing with the spatial coordinate:
  that is, increasing from left to right (if \code{vertical=FALSE})
  or increasing from bottom to top (if \code{vertical=TRUE}).
  If \code{increasing=FALSE}, this ordering will be reversed.
  The default is \code{increasing=TRUE} in all cases
  except when \code{vertical=TRUE} and
  the domain of the colourmap is a finite set of discrete inputs.
}
\value{
  None.
}
\seealso{\code{\link{colourmap}}}
\examples{
  co <- colourmap(rainbow(100), breaks=seq(-1,1,length=101))
  plot(co)
  plot(co, col.ticks="pink")
  ca <- colourmap(rainbow(8), inputs=letters[1:8])
  plot(ca, vertical=TRUE)
  ## logarithmic colour map
  cl <- colourmap(rainbow(25), range=c(0.1, 1000), compress=log10)
  plot(cl)
}
\author{
  \adrian
  and \rolf
}
\keyword{spatial}
\keyword{color}
\keyword{hplot}
