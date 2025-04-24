\name{plot.symbolmap}
\alias{plot.symbolmap}
\title{
  Plot a Graphics Symbol Map
}
\description{
  Plot a representation of a graphics symbol map,
  similar to a plot legend.
}
\usage{
\method{plot}{symbolmap}(x, \dots, main, xlim = NULL, ylim = NULL,
                         vertical = FALSE,
                         side = c("bottom", "left", "top", "right"),
                         annotate = TRUE, labelmap = NULL, add = FALSE,
                         nsymbols = NULL, warn = TRUE,
                         colour.only=FALSE,
                         representatives=NULL)
} 
\arguments{
  \item{x}{
    Graphics symbol map (object of class \code{"symbolmap"}).
  }
  \item{\dots}{
    Additional graphics arguments passed to
    \code{\link{points}}, \code{\link{symbols}} or \code{\link{axis}}.
  }
  \item{main}{
    Main title for the plot. A character string.
  }
  \item{xlim,ylim}{
    Coordinate limits for the plot.
    Numeric vectors of length 2.
  }
  \item{vertical}{
    Logical. Whether to plot the symbol map in a vertical orientation.
  }
  \item{side}{
    Character string specifying the position of the text
    that annotates the symbols. Alternatively one of the integers 1 to 4.
  }
  \item{annotate}{
    Logical. Whether to annotate the symbols with labels.
  }
  \item{labelmap}{
    Transformation of the labels.
    A function or a scale factor which will be applied to
    the data values corresponding to the plotted symbols.
  }
  \item{add}{
    Logical value indicating whether to add the plot to the
    current plot (\code{add=TRUE}) or to initialise a new plot.
  }
  \item{nsymbols}{
    Optional. The maximum number of symbols that should be displayed.
    Ignored if \code{representatives} are given.
  }
  \item{warn}{
    Logical value specifying whether to issue a warning
    when the plotted symbol map does not represent every possible
    discrete value.
  }
  \item{colour.only}{
    Logical value. If \code{TRUE}, the colour map information will be
    extracted from the symbol map, and only this colour map will be
    plotted. If \code{FALSE} (the default) the entire symbol map is
    plotted, including information about symbol shape and size as well
    as colour.
  }
  \item{representatives}{
    Optional. Vector containing the values of the input data
    which should be shown on the plot.
  }
}
\details{
  A graphics symbol map (object of class \code{"symbolmap"})
  is an association between data values and graphical symbols. 
  
  This command plots the graphics symbol map itself,
  in the style of a plot legend.

  For a map of continuous values (a symbol map which represents a range
  of numerical values) the plot will select about \code{nsymbols}
  different values within this range, and plot their graphical
  representations.

  For a map of discrete inputs (a symbol map which represents a finite
  set of elements, such as categorical values) the plot will try to
  display the graphical representation of every possible input,
  up to a maximum of \code{nsymbols} items. If there are more than
  \code{nsymbols} possible inputs, a warning will be issued (if
  \code{warn=TRUE}, the default).
}
\section{Nonlinear transformation}{
  For a map of continuous values, if \code{representatives} is not
  specified, then the default is to choose approximately \code{nsymbols}
  values equally spaced within the numerical range.

  However, if the symbol map includes a nonlinear transformation
  (arguments \code{compress} and \code{decompress} were passed to
  \code{\link{symbolmap}}) or if one of its graphical parameters
  is a colour map which includes a nonlinear transformation 
  (arguments \code{compress} and \code{decompress} were passed to
  \code{\link{colourmap}}) then the representative values will
  be equally spaced on the compressed scale. 
}
\value{
  None.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{symbolmap}} to create a symbol map.
  
  \code{\link{invoke.symbolmap}} to apply the symbol map to some data
  and plot the resulting symbols.
}
\examples{
  g <- symbolmap(inputs=letters[1:10], pch=11:20)
  plot(g)

  ## squares and circles with area proportional to |x|
  ## For 'squares', size is side length; size = sqrt(area)
  ## For 'circles', size is diameter;    size = sqrt(area * 4/pi)
  g2 <- symbolmap(range=c(-1,1),
                  shape=function(x) ifelse(x > 0, "circles", "squares"),
                  size=function(x) ifelse(x > 0,
                                          sqrt(abs(x)*4/pi),
                                          sqrt(abs(x))),
                  bg = function(x) ifelse(abs(x) < 1, "red", "black"))
  plot(g2, vertical=TRUE, side="left", col.axis="blue", cex.axis=2)
  plot(g2, representatives=c(-1,0,1))

  ## logarithmic display scale
  gl <- symbolmap(range=c(1,1000), pch=21,
                  bg=function(x) ifelse(x < 50, "red", "blue"),
                  cex=function(x) ifelse(x < 100, 1, 2),
                  compress=log10)
  gl(10)
  plot(gl, nsymbols=4)

  ## logarithmic symbol map
  gll <- symbolmap(range=c(1,1000),
                   transform=log10, size=function(x) { 1 + x },
                   compress=log10)
  gll(10)
  plot(gll, nsymbols=4)
}
\keyword{spatial}
\keyword{hplot}
