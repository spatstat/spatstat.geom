\name{plot.ppp}
\alias{plot.ppp}
\title{plot a Spatial Point Pattern}
\description{
  Plot a two-dimensional spatial point pattern
}
\usage{
 \method{plot}{ppp}(x, main, \dots, clipwin=NULL,
                    chars=NULL, cols=NULL,
                    use.marks=TRUE, which.marks=NULL,
                    add=FALSE, type=c("p","n"),
                    legend=TRUE,
                    leg.side=c("left", "bottom", "top", "right"),
                    leg.args=list(),
                    symap=NULL, maxsize=NULL, meansize=NULL, markscale=NULL,
                    minsize=NULL, zerosize=NULL,
                    zap=0.01,
                    show.window=show.all, show.all=!add, do.plot=TRUE,
                    multiplot=TRUE,
                    background=NULL, clip.background=FALSE,
                    scramble.cols=FALSE)
}
\arguments{
  \item{x}{
    The spatial point pattern to be plotted.
    An object of class \code{"ppp"},
    or data which can be converted into 
    this format by \code{\link{as.ppp}()}.
  }
  \item{main}{
    text to be displayed as a title above the plot.
  }
  \item{\dots}{
    extra arguments that will be passed to the plotting functions
    \code{\link{plot.default}}, \code{\link{points}} and/or
    \code{\link{symbols}}. Not all arguments will be recognised.
  }
  \item{clipwin}{
    Optional. A window (object of class \code{"owin"}).
    Only this subset of the image will be displayed.
  }
  \item{chars}{
    the plotting character(s) used to plot points.
    Either a single character, an integer,
    or a vector of single characters or integers.
    Ignored if \code{symap} is given.
  }
  \item{cols}{
    the colour(s) used to plot points.
    Either an integer index from 1 to 8 (indexing the standard colour
    palette), a character string giving the name of a colour,
    or a string giving the hexadecimal representation of a colour,
    or a vector of such integers or strings.
    See the section on \emph{Colour Specification}
    in the help for \code{\link[graphics]{par}}.
    Ignored if \code{symap} is given.
  }
  \item{use.marks}{
    logical flag; if \code{TRUE}, plot points using a different
    plotting symbol for each mark;
    if \code{FALSE}, only the locations of the points will be plotted,
    using \code{\link{points}()}.
  }
  \item{which.marks}{
    Index determining which column of marks to use,
    if the marks of \code{x} are a data frame.
    A character or integer vector identifying one or more
    columns of marks. 
    If \code{add=FALSE} then
    the default is to plot all columns of marks, in a series of separate
    plots.
    If \code{add=TRUE} then only one column of marks can be plotted,
    and the default is \code{which.marks=1}
    indicating the first column of marks.
  }
  \item{add}{
    logical flag; if \code{TRUE},
    just the points are plotted, over the existing plot.
    A new plot is not created, and
    the window is not plotted.
  }
  \item{type}{
    Type of plot: either \code{"p"} or \code{"n"}.
    If \code{type="p"} (the default), both the points and the observation window
    are plotted. If \code{type="n"}, only the window is plotted.
  }
  \item{legend}{
    Logical value indicating whether to add a legend showing the mapping
    between mark values and graphical symbols (for a marked point pattern).
  }
  \item{leg.side}{
    Position of legend relative to main plot.
  }
  \item{leg.args}{
    List of additional arguments passed to \code{\link{plot.symbolmap}}
    or \code{\link{symbolmap}} to control the legend.
    In addition to arguments documented under
    \code{\link{plot.symbolmap}}, and graphical arguments recognised
    by \code{\link{symbolmap}}, the list may also include the argument
    \code{sep} giving the separation between the main plot and the
    legend, or \code{sep.frac} giving the separation as a fraction
    of the largest dimension (maximum of width and height) of the main
    plot.
  }
  \item{symap}{
    The graphical symbol map to be applied to the marks.
    An object of class \code{"symbolmap"}; see
    \code{\link{symbolmap}}.
  }
  \item{maxsize}{
    \emph{Maximum} physical size of the circles/squares plotted
    when \code{x} is a marked point pattern with 
    numerical marks.
    Incompatible with \code{meansize} and \code{markscale}.
    Ignored if \code{symap} is given.
  }
  \item{meansize}{
    \emph{Average} physical size of the circles/squares plotted
    when \code{x} is a marked point pattern with 
    numerical marks.
    Incompatible with \code{maxsize} and \code{markscale}.
    Ignored if \code{symap} is given.
  }
  \item{markscale}{
    physical scale factor determining the sizes of the
    circles/squares plotted when \code{x} is a marked point pattern with 
    numerical marks. Mark value will be multiplied by \code{markscale}
    to determine physical size.
    Incompatible with \code{maxsize} and \code{meansize}.
    Ignored if \code{symap} is given.
  }
  \item{minsize}{
    \emph{Minimum} physical size of the circles/squares plotted
    when \code{x} is a marked point pattern with 
    numerical marks.
    Incompatible with \code{zerosize}.
    Ignored if \code{symap} is given.
  }
  \item{zerosize}{
    Physical size of the circle/square representing a mark value of zero,
    when \code{x} is a marked point pattern with numerical marks.
    Incompatible with \code{minsize}.
    Defaults to zero.
    Ignored if \code{symap} is given.
  }
  \item{zap}{
    Fraction between 0 and 1.
    When \code{x} is a marked point pattern with numerical marks,
    \code{zap} is the smallest mark value
    (expressed as a fraction of the maximum possible mark) that will
    be plotted. 
    Any points which have marks smaller in absolute value
    than \code{zap * max(abs(marks(x)))} will not be plotted.
  }
  \item{show.window}{
    Logical value indicating whether to plot the observation
    window of \code{x}.
  }
  \item{show.all}{
    Logical value indicating whether to plot everything
    including the main title and the observation window of \code{x}.
  }
  \item{do.plot}{
    Logical value determining whether to actually perform the plotting.
  }
  \item{multiplot}{
    Logical value giving permission to display multiple plots.
  }
  \item{background}{
    Optional. An object to be plotted underneath the main plot.
    See the section on \emph{Background}.
  }
  \item{clip.background}{
    Logical value specifying whether \code{background} should be
    intersected with the rectangular frame of \code{x}.
  }
  \item{scramble.cols}{
    Logical value. If \code{TRUE}, the sequence of colour values
    will be randomly permuted. This is a useful trick when it is
    desired that adjacent colours in the sequence should be easily
    distinguishable. Applicable only when \code{x} is a multitype
    point pattern.
  }
}
\value{
  (Invisible) object of class \code{"symbolmap"}
  giving the correspondence between 
  mark values and plotting characters.
}
\details{
  This is the \code{plot} method for 
  point pattern datasets (of class \code{"ppp"}, see \code{\link{ppp.object}}).

  First the observation window \code{Window(x)} is plotted
  (if \code{show.window=TRUE}).
  Then the points themselves are plotted,
  in a fashion that depends on their marks,
  as follows.
  \describe{
    \item{unmarked point pattern:}{
      If the point pattern does not have marks, or if \code{use.marks = FALSE},
      then the locations of all points will be plotted 
      using a single plot character
    }
    \item{multitype point pattern:}{
      If \code{marks(x)} is a factor, then 
      each level of the factor is 
      represented by a different plot character.
    }
    \item{continuous marks:}{
      If \code{marks(x)} is a numeric vector,
      the marks are rescaled to the unit interval and
      each point is represented by a circle 
      with \emph{diameter} proportional to the rescaled mark
      (if the value is positive) or a square with \emph{side length}
      proportional to the absolute value of the rescaled mark
      (if the value is negative).
    }
    \item{other kinds of marks:}{
      If \code{marks(x)} is neither numeric nor a factor,
      then each possible mark will be represented by a
      different plotting character. The default is to
      represent the \eqn{i}th smallest mark value by
      \code{points(..., pch=i)}.
    }
  }

  If there are several columns of marks, and if \code{which.marks} is
  missing or \code{NULL}, then
  \itemize{
    \item 
    if \code{add=FALSE} and \code{multiplot=TRUE}
    the default is to plot all columns of marks, in a series of separate
    plots, placed side-by-side. The plotting is coordinated by
    \code{\link{plot.listof}}, which calls \code{plot.ppp} to make each of
    the individual plots.
    \item
    Otherwise, only one column of marks can be plotted,
    and the default is \code{which.marks=1}
    indicating the first column of marks.
  }

  Plotting of the window \code{Window(x)} is performed by
  \code{\link{plot.owin}}. This plot may be modified 
  through the \code{...} arguments. In particular the
  extra argument \code{border} determines
  the colour of the window, if the window is not a binary mask.

  Plotting of the points themselves is performed
  by the function \code{\link{points}}, except for the case of
  continuous marks, where it is performed by \code{\link{symbols}}.
  Their plotting behaviour may be modified through the \code{...}
  arguments.

  If the argument \code{symap} is given, then it determines the
  graphical display of the points. It should be a symbol map
  (object of class \code{"symbolmap"}) created by the function
  \code{\link{symbolmap}}.

  If \code{symap} is not given, then the
  following arguments can be used to specify how the points are plotted:
  \itemize{
    \item
    The argument \code{chars} determines the plotting character
    or characters used to display the points (in all cases except
    for the case of continuous marks). For an unmarked point pattern,
    this should be a single integer or character determining a
    plotting character (see \code{par("pch")}).
    For a multitype point pattern, \code{chars} should be a vector
    of integers or characters, of the same length
    as \code{levels(marks(x))}, and then the \eqn{i}th level or type 
    will be plotted using character \code{chars[i]}.
    \item
    If \code{chars} is absent, but there is an extra argument
    \code{pch}, then this will determine the plotting character for
    all points.
    \item
    The argument \code{cols} determines the colour or colours used to
    display the points. For an unmarked point pattern,
    \code{cols} should be a character string
    determining a colour. For a multitype point pattern, \code{cols}
    should be a character vector, of the same length
    as \code{levels(marks(x))}: that is, there is one colour for each
    possible mark value. The \eqn{i}th level or type will
    be plotted using colour \code{cols[i]}. For a point pattern with
    continuous marks, \code{cols} can be either a character string
    or a character vector specifying colour values: the range of mark
    values will be mapped to the specified colours.
    \item
    If \code{cols} is absent, the colours used to plot the
    points may be determined by the extra argument \code{fg}
    (for multitype point patterns) or the extra argument \code{col}
    (for all other cases). Note that specifying \code{col} will also
    apply this colour to the window itself.
    \item
    The default colour for the points is a semi-transparent grey,
    if this is supported by the plot device. This behaviour can be
    suppressed (so that the default colour is non-transparent)
    by setting \code{spatstat.options(transparent=FALSE)}.
    \item 
    The arguments \code{maxsize}, \code{meansize} and \code{markscale} 
    are incompatible with each other (and incompatible with
    \code{symap}).
    The arguments \code{minsize} and \code{zerosize} are incompatible
    with each other (and incompatible with \code{symap}).
    Together, these arguments control the physical size of the circles and
    squares which represent the marks in a point pattern with continuous
    marks. The size of a circle is defined as its \emph{diameter};
    the size of a square is its side length.
    If \code{markscale} is given, then a mark value of \code{m}
    is plotted as a circle of diameter \code{m * markscale + zerosize}
    (if \code{m} is positive) or a square of side
    \code{abs(m) * markscale + zerosize}
    (if \code{m} is negative). If \code{maxsize} is given, then the
    largest mark in absolute value, \code{mmax=max(abs(marks(x)))},
    will be scaled to have physical size \code{maxsize}.
    If \code{meansize} is given, then the
    average absolute mark value, \code{mmean=mean(abs(marks(x)))},
    will be scaled to have physical size \code{meansize}.
    If \code{minsize} is given, then the
    minimum mark value, \code{mmean=mean(abs(marks(x)))},
    will be scaled to have physical size \code{minsize}.
    \item
    The user can set the default values of these plotting parameters
    using \code{\link{spatstat.options}("par.points")}.
  }

  To zoom in (to view only a subset of the point pattern at higher
  magnification), use the graphical arguments
  \code{xlim} and \code{ylim} to specify the rectangular field of view.

  The value returned by this plot function is an object of
  class \code{"symbolmap"} representing the mapping from mark values
  to graphical symbols. See \code{\link{symbolmap}}.
  It can be used to make a suitable legend,
  or to ensure that two plots use the same graphics map.
}
\section{Layout of the plot}{
  \itemize{
    \item
    \bold{Placement of main title:}
    The left-right placement of the main title is controlled by the
    argument \code{adj.main} passed to \code{\link[spatstat.geom]{plot.owin}}.
    To remove the main title, set \code{main=""}.
    \item
    \bold{Removing White Space Around The Plot:}
    A frequently-asked question is: How do I remove the white space around
    the plot? Currently \code{plot.ppp} uses the base graphics system of
    \R, so the space around the plot is controlled by parameters 
    to \code{\link{par}}. To reduce the white space, change the
    parameter \code{mar}. Typically, \code{par(mar=rep(0.5, 4))} is
    adequate, if there are no annotations or titles outside the window.
    \item \bold{Drawing coordinate axes and axis labels:}
    Coordinate axes and axis labels are not drawn, by default.
    To draw coordinate axes, set \code{axes=TRUE}.
    To draw axis labels, set \code{ann=TRUE} and give values to the
    arguments \code{xlab} and \code{ylab}. See the Examples.
    Only the default style of axis is supported;
    for more control over the placement and style of axes,
    use the graphics commands
    \code{\link[graphics]{axis}} and \code{\link[graphics]{mtext}}.
  }
}
\section{The Symbol Map}{
  The behaviour of \code{plot.ppp} is different from the
  behaviour of the base \R graphics functions
  \code{\link[graphics]{points}} and \code{\link[graphics]{symbols}}.

  In the base graphics functions \code{\link[graphics]{points}} and \code{\link[graphics]{symbols}},
  arguments such as \code{col}, \code{pch} and \code{cex}
  can be vectors which specify the \emph{representation of each successive point}.
  For example \code{col[3]} would specify the colour of the third point
  in the sequence of points.  If there are 100 points then
  \code{col} should be a vector of length 100.

  In the \pkg{spatstat} function \code{plot.ppp},
  arguments such as \code{col}, \code{pch} and \code{cex} specify
  the \emph{mapping from point characteristics to graphical parameters}
  (called the symbol map).
  For example \code{col[3]} specifies the colour of the third \bold{type of point} in a
  pattern of points of different types. If there are 4 types of points
  then \code{col} should be a vector of length 4.
  
  To modify a symbol map, for example to change the colours used
  without changing anything else, use \code{\link{update.symbolmap}}.
}
\section{Scrambling the colours}{
  When \code{x} is a multitype point pattern and \code{cols} is given,
  the argument \code{scramble.cols=TRUE}
  will cause the sequence of colours in the colour map to be
  randomly permuted. This is a useful trick when it is desired that
  adjacent colours in the sequence should be easily distinguishable.
  The return value of \code{plot.ppp} will include the randomised colour
  sequence rather than the original colour sequence.

  Randomisation implies that the appearance of the plot,
  and the return value, will be different
  each time the command is executed. For reproducible results,
  set the random generator seed using \code{\link{set.seed}}.
}
\section{Background}{
  The argument \code{background} specifies an object that will be
  plotted before the point pattern \code{x} is plotted,
  and will therefore appear
  underneath it. The background can be any spatial object in the
  \pkg{spatstat} package family. To specify graphical arguments
  for the background, make it into a \code{\link{layered}} object
  using \code{\link{layered}} or \code{\link{as.layered}},
  and assign graphical parameters to it using \code{\link{layerplotargs}}.

  Alternatively \code{background} can be a single colour value,
  indicating that the plot should have this background colour.
}
\section{Troubleshooting: text is clipped}{
  Text that is intended to be drawn in the margins of the plot
  can sometimes be clipped so that it is only partly visible
  or is completely missing. Clipping is not an artefact of
  \pkg{spatstat}; it is a feature of base graphics in R.
  To prevent clipping of text, first ensure that the plot margins
  are large enough by setting \code{\link{par}('mar')}.
  In base graphics, some plot functions have the side-effect that they
  re-set the clipping region to be a sub-rectangle of the plotting
  region. To override this and ensure that the clipping region
  is as large as possible,
  set \code{par(xpd=NA)} immediately before calling \code{plot.ppp}.
}
\seealso{
  \code{\link{ppp.object}},
  \code{\link{default.symbolmap.ppp}}.
  
  \code{\link{plot}},
  \code{\link{par}},
  \code{\link{points}},
  \code{\link{text.ppp}},
  \code{\link{plot.owin}},
  \code{\link{symbols}}.

  See also the command \code{iplot} in the \pkg{spatstat.gui} package.
}
\examples{
   plot(cells)

   plot(cells, pch=16)

   # make the plotting symbols larger (for publication at reduced scale)
   plot(cells, cex=2)

   # set it in spatstat.options
   oldopt <- spatstat.options(par.points=list(cex=2))
   plot(cells)
   spatstat.options(oldopt)

   # multitype 
   plot(lansing)

   # marked by a real number
   plot(longleaf)

   # just plot the points
   plot(longleaf, use.marks=FALSE)
   plot(unmark(longleaf)) # equivalent

   # point pattern with multiple marks
   plot(finpines)
   plot(finpines, which.marks="height")

   # controlling COLOURS of points
   plot(cells, cols="blue")
   plot(lansing, cols=c("black", "yellow", "green", 
                        "blue","red","pink"))
   plot(longleaf, fg="blue")
   plot(lansing, cols=rainbow, scramble.cols=TRUE)

   # make window purple
   plot(lansing, border="purple")
   # make everything purple
   plot(lansing, border="purple", cols="purple", col.main="purple",
                 leg.args=list(col.axis="purple"))
 
   # controlling PLOT CHARACTERS for multitype pattern
   plot(lansing, chars = 11:16)
   plot(lansing, chars = c("o","h","m",".","o","o"))

   ## multitype pattern mapped to symbols
   plot(amacrine, shape=c("circles", "squares"), size=0.04)
   plot(amacrine, shape="arrows", direction=c(0,90), size=0.07)

   ## plot trees as trees!
   plot(lansing, shape="arrows", direction=90, cols=1:6)

   # controlling MARK SCALE for pattern with numeric marks
   plot(longleaf, markscale=0.1)
   plot(longleaf, maxsize=5)
   plot(longleaf, meansize=2)
   plot(longleaf, minsize=2)

   # draw circles of diameter equal to nearest neighbour distance
   plot(cells \%mark\% nndist(cells), markscale=1, legend=FALSE)

   # inspecting the symbol map
   v <- plot(amacrine)
   v

   ## variable colours ('cols' not 'col')
   plot(longleaf, cols=function(x) ifelse(x < 30, "red", "black"))

   ## re-using the same mark scale
   a <- plot(longleaf)
   juveniles <- longleaf[marks(longleaf) < 30]
   plot(juveniles, symap=a)

   ## numerical marks mapped to symbols of fixed size with variable colour
   ra <- range(marks(longleaf))
   colmap <- colourmap(terrain.colors(20), range=ra)
   ## filled plot characters are the codes 21-25
   ## fill colour is indicated by 'bg'
   ## outline colour is 'fg'
   sy <- symbolmap(pch=21, bg=colmap, fg=colmap, range=ra)
   plot(longleaf, symap=sy)

   ## or more compactly..
   plot(longleaf, bg=terrain.colors(20), pch=21, cex=1)

   ## plot only the colour map (since the symbols have fixed size and shape)
   plot(longleaf, symap=sy, leg.args=list(colour.only=TRUE))

   ## plot with a background
   plot(bei \%mark\% nndist(bei), background=bei.extra$elev)

   ## clipping
   plot(humberside)
   B <- owin(c(4810, 5190), c(4180, 4430))
   plot(B, add=TRUE, border="red")
   plot(humberside, clipwin=B, main="Humberside (clipped)")

   ## coordinate axes and labels
   plot(humberside, axes=TRUE)
   plot(humberside,            ann=TRUE, xlab="Easting", ylab="Northing")
   plot(humberside, axes=TRUE, ann=TRUE, xlab="Easting", ylab="Northing")
}

\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{hplot}
