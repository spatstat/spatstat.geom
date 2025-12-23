\name{connected}  
\Rdversion{1.1}
\alias{connected}
\alias{connected.im}
\alias{connected.owin}
\title{
  Connected components
}
\description{
  Finds the topologically-connected components of a spatial object,
  such as the connected clumps of pixels in a binary image.
}
\usage{
connected(X, \dots)

\method{connected}{owin}(X, \dots, polygonal=FALSE, method="C", connect=8)

\method{connected}{im}(X, \dots, background, method="C", connect=8)
}
\arguments{
  \item{X}{
    A spatial object such as a pixel image (object of class \code{"im"})
    or a window (object of class \code{"owin"}).
  }
  \item{polygonal}{
    Logical value specifying whether to use polygonal geometry.
    See Details.
  }
  \item{background}{
    Optional. Treat pixels with this value 
    as being part of the background.
  }
  \item{method}{
    String indicating the algorithm to be used. Either \code{"C"}
    or \code{"interpreted"}. See Details.
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat.geom]{as.mask}} to determine the
    pixel resolution.
  }
  \item{connect}{
    The connectivity of the pixel grid: either 8 or 4.
  }
}
\details{
  The function \code{connected} is generic, with methods
  for pixel images (class \code{"im"}) and windows (class \code{"owin"})
  described here. There are also methods for
  tessellations (\code{\link{connected.tess}}),
  point patterns (\code{\link{connected.ppp}} and
  \code{connected.lpp}),
  and linear networks (\code{connected.linnet}).
  
  The functions described here compute the connected component transform
  (Rosenfeld and Pfalz, 1966)
  of a binary image or binary mask.
  The pixel values in \code{X} are first mapped to logical values.
  Then the algorithm
  identifies the connected components (topologically-connected clumps
  of pixels) in the foreground.

  When \code{X} is converted to a logical valued image,
  the pixel value \code{NA} is always mapped to \code{FALSE}.
  If \code{background} is specified, then the value \code{background} is also
  mapped to \code{FALSE}.   All other values are mapped to \code{TRUE}.
  In the special case where \code{X} is already a logical-valued
  image, the default is \code{background=FALSE}, so that the
  values \code{NA} and \code{FALSE} are mapped to \code{FALSE}.

  Two pixels belong to the same connected component if they have the value
  \code{TRUE} and if they are neighbours.
  This rule is applied repeatedly until it terminates.
  Then each connected component
  contains all the pixels that can be reached by stepping from neighbour
  to neighbour.

  Pixels are defined to be neighbours if they are physically adjacent to
  each other. If \code{connect=4}, each pixel has 4 neighbours,
  lying one step above or below, or one step to the left or right.  
  If \code{connect=8} (the default), each pixel has 8 neighbours,
  lying one step above or below, or one step to the left or right,
  or one diagonal step away. (Pixels at the edge of the image have fewer
  neighbours.) The 8-connected algorithm is the default
  because it gives better results when the pixel grid is coarse. The
  4-connected algorithm is faster and is recommended when the pixel grid
  is fine.
  
  If \code{method="C"}, the computation is performed by a compiled C language
  implementation of the classical algorithm of Rosenfeld and Pfalz
  (1966). If \code{method="interpreted"}, the computation is performed
  by an \R implementation of the algorithm of Park et al (2000). 

  By default, the result is a factor-valued image,
  with levels that correspond to
  the connected components. The Examples show how to extract each
  connected component as a separate window object.

  If \code{X} is a window and \code{polygonal=TRUE},
  the result is a tessellation
  (object of class \code{"tess"}) whose tiles are the connected components.
}
\value{
  A pixel image (object of class \code{"im"}) with factor values.
  The levels of the factor correspond to the connected components.

  For \code{connected.owin}, if \code{polygonal=TRUE},
  the result is a tessellation
  (object of class \code{"tess"}) whose tiles are the connected components.
}
\references{
  Park, J.-M., Looney, C.G. and Chen, H.-C. (2000)
  Fast connected component labeling algorithm using a divide and conquer
  technique. Pages 373-376 in
  S.Y. Shin (ed) \emph{Computers and Their Applications:} Proceedings of
  the ISCA 15th International Conference on Computers and Their
  Applications, March 29-31, 2000, New Orleans, Louisiana USA. ISCA
  2000, ISBN 1-880843-32-3. 

  Rosenfeld, A. and Pfalz, J.L. (1966)
  Sequential operations in digital processing.
  \emph{Journal of the Association for Computing Machinery} \bold{13}
  471-494.
}
\seealso{
  \code{\link{connected.ppp}},
  \code{\link{connected.tess}},
  \code{\link{im.object}}, 
  \code{\link{tess}}
}
\section{Warnings}{
  It may be hard to distinguish different components 
  in the default plot because the colours of nearby components may be
  very similar. See the Examples for a randomised colour map.
  
  The algorithm for \code{method="interpreted"}
  can be very slow for large images (or images where
  the connected components include a large number of pixels).
}
\examples{
  d <- distmap(cells, dimyx=256)
  X <- levelset(d, 0.07)
  plot(X)
  Z <- connected(X)
  plot(Z)
  # or equivalently
  Z <- connected(d <= 0.07)

  # number of components
  nc <- length(levels(Z))
  # plot with randomised colour map
  plot(Z, col=hsv(h=sample(seq(0,1,length=nc), nc)))

  # how to extract the components as a list of windows
  W <- tiles(tess(image=Z))

  ## polygonal algorithm
  A <- regularpolygon(7)
  B <- regularpolygon(13, 0.3)
  D <- owin(c(-1, -0.6), c(1, 1.4))
  E <- regularpolygon(5, 0.3)
  Funky <- union.owin(D, E, setminus.owin(A, B))
  plot(Funky, col="blue")
  plot(connected(Funky, polygonal=TRUE), do.col=TRUE)
}
\author{
  Original \R code by Julian Burgos, University of Washington.
  Adapted for \pkg{spatstat} by
  \adrian
  and \rolf.
}
\keyword{spatial}
\keyword{math}
