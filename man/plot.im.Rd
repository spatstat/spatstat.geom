\name{plot.im}
\alias{plot.im}
\alias{image.im}
\title{Plot a Pixel Image}
\description{
  Plot a pixel image.
}
\usage{
  \method{plot}{im}(x, \dots,
                   main,
                   add=FALSE, clipwin=NULL,
                   col=NULL, reverse.col=FALSE,
                   valuesAreColours=NULL, log=FALSE,
                   ncolours=256, gamma=1, 
                   ribbon=show.all, show.all=!add,
                   drop.ribbon=FALSE,
                   ribside=c("right", "left", "bottom", "top"),
                   ribsep=0.15, ribwid=0.05, ribn=1024,
                   ribscale=1, ribargs=list(), riblab=NULL, colargs=list(),
                   useRaster=NULL, workaround=FALSE, zap=1,
                   do.plot=TRUE,
                   addcontour=FALSE, contourargs=list(),
                   background=NULL, clip.background=FALSE)

  \method{image}{im}(x, \dots,
                   main,
                   add=FALSE, clipwin=NULL,
                   col=NULL, reverse.col=FALSE,
                   valuesAreColours=NULL, log=FALSE, 
                   ncolours=256, gamma=1, 
                   ribbon=show.all, show.all=!add,
                   drop.ribbon=FALSE,
                   ribside=c("right", "left", "bottom", "top"),
                   ribsep=0.15, ribwid=0.05, ribn=1024,
                   ribscale=1, ribargs=list(), riblab=NULL, colargs=list(),
                   useRaster=NULL, workaround=FALSE, zap=1,
                   do.plot=TRUE,
                   addcontour=FALSE, contourargs=list(),
                   background=NULL, clip.background=FALSE)
}
\arguments{
  \item{x}{
    The pixel image to be plotted.
    An object of class \code{"im"} (see \code{\link{im.object}}).
  }
  \item{\dots}{
    Extra arguments passed to \code{\link[graphics]{image.default}}
    to control the plot. See Details.
  }
  \item{main}{Main title for the plot.}
  \item{add}{
    Logical value indicating whether to superimpose the image on the
    existing plot (\code{add=TRUE})
    or to initialise a new plot (\code{add=FALSE}, the default).
  }
  \item{clipwin}{
    Optional. A window (object of class \code{"owin"}).
    Only this subset of the image will be displayed.
  }
  \item{col}{
    Colours for displaying the pixel values.
    Either a character vector of colour values,
    an object of class \code{\link{colourmap}},
    or a \code{function} as described under Details.
  }
  \item{reverse.col}{
    Logical value. If \code{TRUE}, the sequence of colour values
    specified by \code{col} will be reversed.
  }
  \item{valuesAreColours}{
    Logical value. If \code{TRUE}, the pixel values of \code{x}
    are to be interpreted as colour values.
  }
  \item{log}{
    Logical value. If \code{TRUE}, the colour map will be
    evenly-spaced on a logarithmic scale.
  }
  \item{ncolours}{
    Integer. The default number of colours in the colour map
    for a real-valued image.
  }
  \item{gamma}{
    Exponent for the gamma correction of the colours.
    A single positive number.
  }
  \item{ribbon}{
    Logical flag indicating whether to display a ribbon
    showing the colour map. Default is \code{TRUE}
    for new plots and \code{FALSE} for added plots.
  }
  \item{show.all}{
    Logical value indicating whether to display all plot elements
    including the main title and colour ribbon. Default is \code{TRUE}
    for new plots and \code{FALSE} for added plots.
  }
  \item{drop.ribbon}{
    Logical value. If \code{TRUE}, then a ribbon will not be displayed
    if all pixel values are equal.
  }
  \item{ribside}{
    Character string indicating where to display the ribbon
    relative to the main image.
  }
  \item{ribsep}{
    Factor controlling the space between the ribbon and the image.
  }
  \item{ribwid}{
    Factor controlling the width of the ribbon.
  }
  \item{ribn}{
    Number of different values to display in the ribbon.
  }
  \item{ribscale}{
    Rescaling factor for tick marks. The values on the numerical scale
    printed beside the ribbon will be multiplied by this rescaling factor.
  }
  \item{ribargs}{
    List of additional arguments passed to
    \code{\link[graphics]{image.default}}, 
    \code{\link[graphics]{axis}} and
    \code{\link[grDevices]{axisTicks}}
    to control the display of the ribbon and its scale axis. These may override
    the \code{\dots} arguments. 
  }
  \item{riblab}{
    Text to be plotted in the margin near the ribbon.
    A character string or expression to be interpreted as text,
    or a list of arguments to be passed to
    \code{\link[graphics]{mtext}}.
  }
  \item{colargs}{
    List of additional arguments passed to
    \code{col} if it is a function.
  }
  \item{useRaster}{
    Logical value, passed to \code{\link[graphics]{image.default}}.
    Images are plotted using a bitmap raster if \code{useRaster=TRUE}
    or by drawing polygons if \code{useRaster=FALSE}.
    Bitmap raster display tends to produce better results,
    but is not supported on all graphics devices.
    The default is to use bitmap raster display if it is supported.
  }
  \item{workaround}{
    Logical value, specifying whether to use a workaround to avoid a bug
    which occurs with some device drivers in \R, in which the image
    has the wrong spatial orientation. See the section on
    \bold{Image is Displayed in Wrong Spatial Orientation} below.
  }
  \item{zap}{
    Noise threshold factor. A numerical value greater than or equal to 1.
    If the range of pixel values is less than
    \code{zap * .Machine$double.eps}, the image will be treated as
    constant. This avoids displaying images which should be constant
    but contain small numerical errors.
  }
  \item{do.plot}{
    Logical value indicating whether to actually plot the image
    and colour ribbon.
    Setting \code{do.plot=FALSE} will simply return the
    colour map and the bounding box that were chosen for the plot.
  }
  \item{addcontour}{
    Logical value specifying whether to add contour lines to the image
    plot. The contour lines will also be drawn on the colour ribbon.
  }
  \item{contourargs}{
    Optional list of arguments to be passed to
    \code{\link[graphics]{contour.default}} to control the contour plot.
  }
  \item{background}{
    Optional. An object to be plotted underneath the image.
    See the section on \emph{Background}.
  }
  \item{clip.background}{
    Logical value specifying whether \code{background} should be
    intersected with the rectangular frame of \code{x}.
  }
}
\value{
  The colour map used. An object of class \code{"colourmap"}.

  Also has an attribute \code{"bbox"} giving a bounding box
  for the plot (containing the main colour image, the
  colour ribbon if plotted, and the background if given).
  If a ribbon was plotted, there is also an attribute
  \code{"bbox.legend"} giving a bounding box for the ribbon image.
  Text annotation occurs outside these bounding boxes.

  Also has an attribute \code{"at"} giving the axis positions
  of the tick marks next to the colour ribbon.
}
\details{
  This is the \code{plot} method for the class \code{"im"}.
  [It is also the \code{image} method for \code{"im"}.]

  The pixel image \code{x} is displayed on the current plot device,
  using equal scales on the \code{x} and \code{y} axes.

  If \code{ribbon=TRUE}, a legend will be plotted.
  The legend consists of a colour ribbon and an axis with tick-marks,
  showing the correspondence between the pixel values and the colour map.

  Arguments \code{ribside, ribsep, ribwid} control the placement
  of the colour ribbon.
  By default, the ribbon is placed at the right of the main image.
  This can be changed using the argument \code{ribside}.
  The width of the ribbon is \code{ribwid} times the size of the pixel
  image, where `size' means the larger of the width and the height.
  The distance separating the ribbon and the image is \code{ribsep} times
  the size of the pixel image.

  The ribbon contains the colours representing \code{ribn}
  different numerical values, evenly spaced between the minimum and
  maximum pixel values in the image \code{x}, rendered according to
  the chosen colour map.

  The argument \code{ribargs} controls the annotation of the
  colour ribbon. It is a list of arguments to be passed to
  \code{\link[graphics]{image.default}}, 
  \code{\link[graphics]{axis}} and
  \code{\link[grDevices]{axisTicks}}.
  To plot the colour ribbon without the axis and
  tick-marks, use \code{ribargs=list(axes=FALSE)}.
  To ensure that the numerals or symbols printed next to the colour map
  are oriented horizontally, use \code{ribargs=list(las=1)}.
  To double the size of the numerals or symbols, 
  use \code{ribargs=list(cex.axis=2)}.
  To control the number of tick-marks, use \code{ribargs=list(nint=N)}
  where \code{N} is the desired number of intervals (so there will
  be \code{N+1} tickmarks, subject to the vagaries of \R internal code). 

  The argument \code{riblab} contains text that should be displayed 
  in the margin next to the ribbon. Alternatively \code{riblab} may be a
  list of arguments to be passed to
  \code{\link[graphics]{mtext}}; in that case \code{riblab$text} should
  specify the text to be displayed, while other arguments
  can specify a different
  location for the text (e.g. \code{riblab$side = 1} for text at the
  bottom of the ribbon) and a different font and colour if desired.
  
  The argument \code{ribscale} is used
  to rescale the numerical values printed next to the colour map,
  for convenience. 
  For example if the pixel values in \code{x} range between 1000 and
  4000, it would be sensible to use \code{ribscale=1/1000} so that the
  colour map tickmarks would be labelled 1 to 4. 

  Normally the pixel values are displayed using the colours given in the
  argument \code{col}. This may be either
  \itemize{
    \item an explicit colour map (an object of class
    \code{"colourmap"}, created by the command \code{\link{colourmap}}).
    This is the best way to ensure
    that when we plot different images, the colour maps are consistent.
    \item a character vector or integer vector
    that specifies a set of colours.
    The colour mapping will be stretched to match the range of
    pixel values in the image \code{x}. The mapping of pixel values
    to colours is determined as follows.
    \describe{
      \item{logical-valued images:}{the values \code{FALSE} and
	\code{TRUE} are mapped to the colours \code{col[1]} and
	\code{col[2]} respectively.
	The vector \code{col} should have length 2.
      }
      \item{factor-valued images:}{the factor levels \code{levels(x)}
	are mapped to the entries of \code{col} in order. The vector
	\code{col} should have the same length as \code{levels(x)}.
      }
      \item{numeric-valued images:}{
	By default, the range of pixel values in \code{x}
	is divided into \code{n = length(col)}
	equal subintervals, which are mapped to the colours in \code{col}.
	(If \code{col} was not specified, it defaults to a vector of 255
	colours.)
      
	Alternatively if the argument \code{zlim} is given, it should be
	a vector of length 2 specifying an interval of real numbers.
	This interval will be used instead of the range of pixel
	values. The interval from \code{zlim[1]} to \code{zlim[2]} will be
	mapped to the colours in \code{col}. This facility enables the user to
	plot several images using a consistent colour map.

	Alternatively if the argument \code{breaks} is given,
	then this specifies
	the endpoints of the subintervals that are mapped to each colour.
	This is incompatible with \code{zlim}.
      
	The arguments \code{col} and \code{zlim} or \code{breaks}
	are then passed to the function \code{\link{image.default}}.
	For examples of the use of these arguments,
	see \code{\link{image.default}}.
      }
    }
    \item {
      a \code{function} in the \R language
      with an argument named \code{range}
      or \code{inputs}.
    
      If \code{col} is a function with an argument named \code{range},
      and if the pixel values of \code{x} are numeric values, 
      then the colour values will be determined by evaluating
      \code{col(range=range(x))}. The result of this evaluation
      should be a character vector containing colour values, or
      a \code{"colourmap"} object. Examples of such functions
      are \code{\link{beachcolours}} and \code{\link{beachcolourmap}}.
  
      If \code{col} is a function with an argument named \code{inputs},
      and if the pixel values of \code{x} are discrete values (integer,
      logical, factor or character), 
      then the colour values will be determined by evaluating
      \code{col(inputs=p)} where \code{p} is the set of possible pixel
      values. The result should be a character vector 
      containing colour values, or a \code{"colourmap"} object.
    }
    \item{
      a \code{function} in the \R language with first argument
      named \code{n}.
      The colour values will be determined by evaluating
      \code{col(n)} where \code{n} is the
      number of distinct pixel values, up to a maximum of 128.
      The result of this evaluation
      should be a character vector containing color values.
      Examples of such functions are
      \code{\link[grDevices]{heat.colors}},
      \code{\link[grDevices]{terrain.colors}},
      \code{\link[grDevices]{topo.colors}} and
      \code{\link[grDevices]{cm.colors}}. 
   }
  }

  If \code{col} is missing or \code{col=NULL}, the default colour values
  are the linear, perceptually uniform colour sequence given by
  \code{\link[spatstat.data]{Kovesi}[[29,"values"]]}.
  
  If \code{spatstat.options("monochrome")} has been set to \code{TRUE}
  then \bold{all colours will be converted to grey scale values}.

  If \code{reverse.col=TRUE}, the sequence of colour values specified
  by \code{col} will be reversed (unless \code{valuesAreColours=TRUE}).
  
  Other graphical parameters controlling the display of both the pixel image
  and the ribbon can be passed through the \code{...} arguments
  to the function \code{\link[graphics]{image.default}}.
  A parameter is handled only if it is one of the following:
  \itemize{
    \item a formal argument of \code{\link[graphics]{image.default}}
    that is operative when \code{add=TRUE}.
    \item one of the parameters
    \code{"main", "asp", "sub", "axes", "xlab", "ylab"}
    described in \code{\link[graphics]{plot.default}}.
    \item one of the parameters
    \code{"ann", "cex", "font", "cex.axis", "cex.lab", "cex.main", "cex.sub",
      "col.axis", "col.lab", "col.main", "col.sub",
      "font.axis", "font.lab", "font.main", "font.sub"}
    described in \code{\link[graphics]{par}}.
    \item the argument \code{box}, a logical value specifying whether
    a box should be drawn.
  }
  Images are plotted using a bitmap raster if \code{useRaster=TRUE}
  or by drawing polygons if \code{useRaster=FALSE}.
  Bitmap raster display (performed by \code{\link[graphics]{rasterImage}})
  tends to produce better results, but is not supported on all graphics devices.
  The default is to use bitmap raster display if it is
  supported according to \code{\link[grDevices]{dev.capabilities}}.

  Alternatively, the pixel values could be directly 
  interpretable as colour values in \R. That is, the pixel values
  could be character strings that represent colours, or
  values of a factor whose levels are character strings representing
  colours.
  \itemize{
    \item If \code{valuesAreColours=TRUE}, then the pixel values will
    be interpreted as colour values and displayed using these colours.
    \item If \code{valuesAreColours=FALSE}, then the pixel values will
    \emph{not} be interpreted as colour values, even if they could be.
    \item If \code{valuesAreColours=NULL}, the algorithm will guess
    what it should do. If the argument \code{col} is given,
    the pixel values will \emph{not} be interpreted as colour values. Otherwise,
    if all the pixel values are strings that represent colours, then
    they will be interpreted and displayed as colours.
  }
  If pixel values are interpreted as colours, 
  the arguments \code{col} and \code{ribbon} will be ignored,
  and a ribbon will not be plotted.
}
\section{Adding contour lines}{
  If \code{addcontour=TRUE}, contour lines will be superimposed on the
  image. Lines will also be superimposed on the colour ribbon at the
  corresponding positions.
  The display of the contour lines can be controlled by
  arguments specified in \code{contourargs}.
}
\section{Main title, box and white space}{
  The left-right placement of the main title is controlled by the
  argument \code{adj.main} passed to \code{\link[spatstat.geom]{plot.owin}}.
  A rectangular box surrounding the image domain is drawn by default; it can be
  suppressed by setting \code{box=FALSE}.
  White space around the plot is controlled by
  \code{\link[graphics]{par}('mar')}.
}
\section{Complex-valued images}{
  If the pixel values in \code{x} are complex numbers,
  they will be converted into four images containing the real and
  imaginary parts and the modulus and argument, 
  and plotted side-by-side using \code{\link{plot.imlist}}.
}
\section{Monochrome colours}{
  If \code{spatstat.options("monochrome")} has been set to \code{TRUE},
  then \bold{the image will be plotted in greyscale}.
  The colours are converted to grey scale values using
  \code{\link{to.grey}}.
  The choice of colour map still has an effect, since it determines
  the final grey scale values.

  Monochrome display can also be achieved by
  setting the graphics device parameter \code{colormodel="grey"}
  when starting a new graphics device, or in a call to
  \code{\link{ps.options}} or \code{\link{pdf.options}}.
}
\section{Background}{
  The argument \code{background} specifies an object that will be
  plotted before the image \code{x} is plotted, and will therefore appear
  underneath it. The background can be any spatial object in the
  \pkg{spatstat} package family. To specify graphical arguments
  for the background, make it into a \code{\link{layered}} object
  using \code{\link{layered}} or \code{\link{as.layered}},
  and assign graphical parameters to it using \code{\link{layerplotargs}}.
}
\section{Troubleshooting: Image Looks Like Noise}{
  An image plot which looks like digital noise
  can be produced when the pixel values are almost exactly equal
  but include a tiny amount of numerical error.
  To check this, look at the numerals plotted next to the colour ribbon,
  or compute \code{diff(range(x))}, to determine whether the range
  of pixel values is almost zero. The behaviour can be suppressed
  by picking a larger value of the argument \code{zap}.
}
\section{Troubleshooting: Image Rendering Errors and Problems}{
  The help for \code{\link[graphics]{image.default}}
  and \code{\link[graphics]{rasterImage}} explains that
  errors may occur, or images may be rendered incorrectly, on some
  devices, depending on the availability of colours and other
  device-specific constraints.

  If the image is not displayed at all,
  try setting \code{useRaster=FALSE} in the call to \code{plot.im}.
  If the ribbon colours are not displayed, set
  \code{ribargs=list(useRaster=FALSE)}.

  Errors may occur on some graphics devices if the image is very
  large. If this happens, try setting \code{useRaster=FALSE} in the
  call to \code{plot.im}.

  The error message
  \code{useRaster=TRUE can only be used with a regular grid}
  means that the \eqn{x} and \eqn{y} coordinates of the pixels in the
  image are not perfectly equally spaced, due to numerical rounding.
  This occurs with some images created by earlier versions of \pkg{spatstat}.
  To repair the coordinates in an image \code{X}, type
  \code{X <- as.im(X)}.
}
\section{Troubleshooting: Image is Displayed in Wrong Spatial Orientation}{
  If the image is displayed in the wrong spatial orientation,
  and you created the image data directly, 
  please check that you understand the \pkg{spatstat} convention for
  the spatial orientation of pixel images. The row index of the matrix
  of pixel values corresponds to the increasing \eqn{y} coordinate; the
  column index of the matrix corresponds to the increasing \eqn{x} coordinate
  (Baddeley, Rubak and Turner, 2015, section 3.6.3, pages 66--67).

  Images can be displayed in the wrong spatial orientation
  on some devices, due to a bug in the device driver. This occurs only
  when the plot coordinates are \emph{reversed}, that is, when 
  the plot was initialised with coordinate limits \code{xlim, ylim}
  such that \code{xlim[1] > xlim[2]} or \code{ylim[1] > ylim[2]} or
  both. This bug is reported to occur only when \code{useRaster=TRUE}.
  To fix this, try setting \code{workaround=TRUE}, or if that is
  unsuccessful, \code{useRaster=FALSE}.
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
  set \code{par(xpd=NA)} immediately before calling \code{plot.im}.
}
\seealso{
  \code{\link{im.object}},
  \code{\link{colourmap}},
  \code{\link{contour.im}},
  \code{\link{persp.im}},
  \code{\link{hist.im}},
  \code{\link[graphics]{image.default}},
  \code{\link{spatstat.options}},
  \code{\link{default.image.colours}}
}
\examples{
   # an image
   Z <- setcov(owin())
   plot(Z)
   plot(Z, ribside="bottom")
   # stretchable colour map
   plot(Z, col=rainbow)
   plot(Z, col=terrain.colors(128), axes=FALSE)
   # fixed colour map
   tc <- colourmap(rainbow(128), breaks=seq(-1,2,length=129))
   plot(Z, col=tc)
   # colour map function, with argument 'range'
   plot(Z, col=beachcolours, colargs=list(sealevel=0.5))
   # tweaking the plot
   plot(Z, main="La vie en bleu", col.main="blue", cex.main=1.5,
        box=FALSE,
        ribargs=list(col.axis="blue", col.ticks="blue", cex.axis=0.75),
        riblab=list(text="couleur", side=1, col="blue", offset=1, font=2))
   # add axes and axis labels
   plot(Z, axes=TRUE, ann=TRUE, xlab="Easting", ylab="Northing")
   # add contour lines
   plot(Z, addcontour=TRUE, contourargs=list(col="white", drawlabels=FALSE))
   # log scale
   V <- eval.im(exp(exp(Z+2))/1e4)
   plot(V, log=TRUE, main="Log scale")
   # it's complex
   Y <- exp(Z + V * 1i)
   plot(Y)
}
\references{
  \baddrubaturnbook
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
 
\keyword{hplot}

