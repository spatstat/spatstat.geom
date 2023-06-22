\name{as.mask}
\alias{as.mask}
\title{Pixel Image Approximation of a Window}
\description{
  Obtain a discrete (pixel image) approximation of a given window
}
\usage{
 as.mask(w, eps=NULL, dimyx=NULL, xy=NULL,
         rule.eps=c("adjust.eps", "grow.frame", "shrink.frame"))
}
\arguments{
  \item{w}{
    A window (object of class \code{"owin"}) or data acceptable
    to \code{\link{as.owin}}.
  }
  \item{eps}{
    (optional) width and height of pixels.
    A single number, or a numeric vector of length 2.
  }
  \item{dimyx}{
    (optional) pixel array dimensions.
    A single integer, or an integer vector of length 2
    giving dimensions in the \emph{y} and \emph{x} directions.
  }
  \item{xy}{
    (optional) data containing pixel coordinates,
    such as a pixel image (object of class \code{"im"}),
    or a window of type \code{"mask"}. See Details.
  }
  \item{rule.eps}{
    Character string (partially matched) specifying what to do when
    \code{eps} is not a divisor of the frame size.
    Ignored if \code{eps} is missing or null. See Details.
  }
}
\value{
  A window (object of class \code{"owin"})
  of type \code{"mask"} representing a binary pixel image.
}
\details{
  A \sQuote{mask} is a spatial window that is
  represented by a pixel image with binary values.
  It is an object of class \code{"owin"} with type \code{"mask"}.
  
  This function \code{as.mask} creates a representation of any spatial window
  \code{w} as a mask. It generates a rectangular grid of
  locations in the plane, tests whether each of these locations lies
  inside \code{w}, and stores the results as a mask.

  The most common use of this function is to approximate the shape
  of a rectangular or polygonal window \code{w} by a mask,
  for computational purposes. In this case,
  we will usually want to have a very fine grid of pixels.

  This function can also be used to generate a coarsely-spaced grid of
  locations inside a window, for purposes such as subsampling
  and prediction.

  The argument \code{w} should be a window (object of class
  \code{"owin"}). If it is another kind of spatial data,
  then the window information will be extracted using \code{\link{as.owin}}. 
  
  The grid spacing and location are controlled by the
  arguments \code{eps}, \code{dimyx} and \code{xy},
  which are mutually incompatible. 

  If \code{eps} is given, then it specifies the \emph{desired} grid spacing,
  that is, the desired size of the pixels. If \code{eps} is a single number,
  it specifies that the desired grid spacing is \code{eps}
  in both the \eqn{x} and \eqn{y} directions, that is, the desired pixels are
  squares with side length \code{eps}. If \code{eps} is a
  vector of length 2, it specifies that the desired grid spacing is
  \code{eps[1]} in the \eqn{x} direction and 
  \code{eps[2]} in the \eqn{y} direction. That is, the desired pixels are
  rectangles of width \code{eps[1]} and height \code{eps[2]}.

  When \code{eps} is given,
  the argument \code{rule.eps} specifies what to do if pixels
  of the desired size would not fit exactly into the
  rectangular frame of \code{w}. 
  \itemize{
    \item
    if \code{rule.eps="adjust.eps"} (the default),
    the rectangular frame will remain unchanged, and
    the grid spacing (pixel size) \code{eps} will be reduced slightly
    so that an integer number of pixels fits exactly into the frame.
    \item if \code{rule.eps="grow.frame"},
    the grid spacing (pixel size) \code{eps} will remain unchanged,
    and the rectangular frame will be expanded slightly
    so that it consists of an integer number of pixels in each direction.
    \item if \code{rule.eps="shrink.frame"},
    the grid spacing (pixel size) \code{eps} will remain unchanged,
    and the rectangular frame will be contracted slightly
    so that it consists of an integer number of pixels in each direction.
  }

  If \code{dimyx} is given, then the pixel grid will be an
  \eqn{m \times n}{m * n} rectangular grid
  where \eqn{m, n} are given by \code{dimyx[2]}, \code{dimyx[1]}
  respectively. \bold{Warning:} \code{dimyx[1]} is the number of
  pixels in the \eqn{y} direction, and \code{dimyx[2]} is the number
  in the \eqn{x} direction. The grid spacing (pixel size) is determined
  by the frame size and the number of pixels.

  If \code{xy} is given, then this should be some kind of
  data specifing the coordinates of a pixel grid. It may be
  \itemize{
    \item 
    a list or structure containing elements \code{x} and \code{y}
    which are numeric vectors of equal length. These will be taken as
    \eqn{x} and \code{y} coordinates of the margins
    of the grid. The pixel coordinates will be generated
    from these two vectors.
    \item
    a pixel image (object of class \code{"im"}).
    \item
    a window (object of class \code{"owin"}) which is
    of type \code{"mask"} so that it contains pixel coordinates.
  }
  If \code{xy} is given and is either a pixel image or a mask,
  then \code{w} may be omitted, and the window
  information will be extracted from \code{xy}.
  
  If neither \code{eps} nor \code{dimyx} nor \code{xy} is given,
  the pixel raster dimensions are obtained from
  \code{\link{spatstat.options}("npixel")}.

  There is no inverse of this function. However, the function
  \code{\link{as.polygonal}} will compute a polygonal approximation
  of a binary mask.
}
\section{Discretisation rule}{
  The rule used in \code{as.mask} is that a pixel is part of the
  discretised window if and only if the centre of the pixel falls in the
  original window. This is usually sufficient for most purposes,
  and is fast to compute.

  Other discretisation rules are possible; they are available using the
  function \code{\link{owin2mask}}.
}
\section{Converting a spatial pattern to a mask}{
  If the intention is to discretise or pixellate a spatial pattern,
  such as a point pattern, line segment pattern or a linear network, then
  \code{as.mask} is not the appropriate function to use, because
  \code{as.mask} extracts only the window information and converts this
  window to a mask.

  To discretise a point pattern,
  use \code{\link{pixellate.ppp}}. To discretise a line segment pattern,
  use \code{\link{pixellate.psp}} or \code{\link{psp2mask}}.
  To discretise a linear network, use
  \code{\link[spatstat.linnet:methods.linnet]{pixellate.linnet}}.
}
\seealso{
  \code{\link{owin2mask}}.
  
  \code{\link{owin.object}},
  \code{\link{as.rectangle}},
  \code{\link{as.polygonal}},
  \code{\link{spatstat.options}}
}
\examples{
  w <- owin(c(0,10),c(0,10), poly=list(x=c(1,2,3,2,1), y=c(2,3,4,6,7)))
  m <- as.mask(w)
  if(interactive()) {
     plot(w)
     plot(m)
  }
  x <- 1:9
  y <- seq(0.25, 9.75, by=0.5)
  m <- as.mask(w, xy=list(x=x, y=y))

  B <- square(1)
  as.mask(B, eps=0.3)
  as.mask(B, eps=0.3, rule.eps="g")
  as.mask(B, eps=0.3, rule.eps="s")
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{manip}
