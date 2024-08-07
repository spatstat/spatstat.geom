\name{as.im}
\alias{as.im}
\alias{as.im.im}
\alias{as.im.owin}
\alias{as.im.matrix}
\alias{as.im.tess}
\alias{as.im.function}
\alias{as.im.funxy}
\alias{as.im.expression}
\alias{as.im.distfun}
\alias{as.im.nnfun}
\alias{as.im.data.frame}
\alias{as.im.default}
\title{Convert to Pixel Image}
\description{
  Converts various kinds of data to a pixel image
}
\usage{
  as.im(X, \dots)

  \method{as.im}{im}(X, W=NULL, \dots,
        eps=NULL, dimyx=NULL, xy=NULL,
        rule.eps=c("adjust.eps", "grow.frame", "shrink.frame"),
        na.replace=NULL)

  \method{as.im}{owin}(X, W=NULL, \dots,
        eps=NULL, dimyx=NULL, xy=NULL,
        rule.eps=c("adjust.eps", "grow.frame", "shrink.frame"),
        na.replace=NULL, value=1)

  \method{as.im}{matrix}(X, W=NULL, \dots)

  \method{as.im}{tess}(X, W=NULL, \dots, 
        eps=NULL, dimyx=NULL, xy=NULL,
        rule.eps=c("adjust.eps", "grow.frame", "shrink.frame"),
        na.replace=NULL, values=NULL)

  \method{as.im}{function}(X, W=NULL, \dots,
        eps=NULL, dimyx=NULL, xy=NULL,
        rule.eps=c("adjust.eps", "grow.frame", "shrink.frame"),
        na.replace=NULL,
        stringsAsFactors=NULL,
        strict=FALSE, drop=TRUE)

  \method{as.im}{funxy}(X, W=Window(X), \dots)

  \method{as.im}{expression}(X, W=NULL, \dots)

  \method{as.im}{distfun}(X, W=NULL, \dots,
        eps=NULL, dimyx=NULL, xy=NULL,
        rule.eps=c("adjust.eps", "grow.frame", "shrink.frame"),
        na.replace=NULL, approx=TRUE)

  \method{as.im}{nnfun}(X, W=NULL, \dots,
        eps=NULL, dimyx=NULL, xy=NULL,
        rule.eps=c("adjust.eps", "grow.frame", "shrink.frame"),
        na.replace=NULL, approx=TRUE)

  \method{as.im}{data.frame}(X, \dots, step, fatal=TRUE, drop=TRUE)

  \method{as.im}{default}(X, W=NULL, \dots,
        eps=NULL, dimyx=NULL, xy=NULL,
        rule.eps=c("adjust.eps", "grow.frame", "shrink.frame"),
        na.replace=NULL)


}
\arguments{
  \item{X}{Data to be converted to a pixel image.}
  \item{W}{Window object which determines the spatial domain
    and pixel array geometry.
  }
  \item{\dots}{Additional arguments passed to \code{X} when \code{X}
    is a function.}
  \item{eps,dimyx,xy,rule.eps}{
    Optional parameters passed to \code{\link[spatstat.geom]{as.mask}}
    which determine the pixel array geometry.
    See \code{\link[spatstat.geom]{as.mask}}.
  }
  \item{na.replace}{Optional value to replace \code{NA} entries in the
    output image.
  }
  \item{value}{Optional.
    The value to be assigned to pixels inside the window,
    if \code{X} is a window. A single atomic value (numeric, integer,
    logical etc).
  }
  \item{values}{Optional.
    Vector of values to be assigned to each tile of the tessellation,
    when \code{X} is a tessellation. An atomic vector (numeric,
    integer, logical etc.)
  }
  \item{strict}{
    Logical value indicating whether to match formal arguments
    of \code{X} when \code{X} is a function. If \code{strict=FALSE}
    (the default), all the \code{\dots} arguments are passed to
    \code{X}. If \code{strict=TRUE}, only named arguments
    are passed, and only if they match the names of formal arguments of
    \code{X}. 
  }
  \item{step}{
    Optional. A single number, or numeric vector of length 2,
    giving the grid step lengths
    in the \eqn{x} and \eqn{y} directions.
  }
  \item{fatal}{
    Logical value indicating what to do if the resulting image
    would be too large for available memory. If \code{fatal=TRUE} (the
    default), an error occurs. If \code{fatal=FALSE}, a warning is
    issued and \code{NULL} is returned.
  }
  \item{drop}{
    Logical value indicating what to do if the result would
    normally be a list of pixel images but the list contains
    only one image.
    If \code{drop=TRUE} (the default), the pixel image is extracted
    and the result is a pixel image.
    If \code{drop=FALSE}, this list is returned as the result.
  }
  \item{stringsAsFactors}{
    Logical value (passed to \code{\link[base]{data.frame}})
    specifying how to handle pixel values which
    are character strings. If \code{TRUE}, character values are
    interpreted as factor levels. If \code{FALSE}, they remain
    as character strings. The default depends on the version of \R.
    See section \emph{Handling Character Strings}.
  }
  \item{approx}{
   Logical value indicating whether to compute an approximate result
   at faster speed.
 }
}
\details{
  This function converts the data \code{X} into a pixel image
  object of class \code{"im"} (see \code{\link{im.object}}).
  The function \code{as.im} is generic, with methods for the classes
  listed above.

  Currently \code{X} may be any of the following:
  \itemize{
    \item
    a pixel image object, of class \code{"im"}.
    \item
    a window object, of class \code{"owin"} (see
    \code{\link{owin.object}}). The result is an image
    with all pixel entries equal to \code{value} inside the window \code{X},
    and \code{NA} outside.
    \item
    a matrix. 
    \item
    a tessellation (object of class \code{"tess"}).
    By default, the result is a factor-valued image, with one factor level
    corresponding to each tile of the tessellation. Pixels are classified
    according to the tile of the tessellation into which they fall.
    If argument \code{values} is given, the result is a pixel image
    in which every pixel inside the \code{i}-th tile of the tessellation
    has pixel value equal to \code{values[i]}.
    \item
    a single number (or a single logical, complex, factor or character
    value). The result is an image
    with all pixel entries equal to this constant value
    inside the window \code{W} (and \code{NA} outside, unless the
    argument \code{na.replace} is given).
    Argument \code{W} is required.
    \item
    a function of the form \code{function(x, y, ...)}
    which is to be evaluated to yield the image pixel values.
    In this case, the additional argument \code{W} must be present.
    This window will be converted to
    a binary image mask. Then the function \code{X} will be evaluated
    in the form \code{X(x, y, ...)} where \code{x} and \code{y} are
    \bold{vectors} containing the \eqn{x} and \eqn{y} coordinates
    of all the pixels in the image mask, and \code{...} are any extra
    arguments given. This function must return a
    vector or factor of the same length as the input vectors,
    giving the pixel values.
    \item
    an object of class \code{"funxy"} representing a
    \code{function(x,y,...)} defined in a spatial region.
    The function will be evaluated as described above.
    The window \code{W} defaults to the domain of definition of the
    function.
    \item
    an object of class \code{"funxy"} which also belongs to one of
    the following special classes. If \code{approx=TRUE} (the default),
    the function will be evaluated approximately using a very fast
    algorithm. If \code{approx=FALSE}, the function will be evaluated
    exactly at each grid location as described above.
    \itemize{
      \item
      an object of class \code{"distfun"} representing a distance function
      (created by the command \code{\link{distfun}}).
      The fast approximation is 
      the distance transform \code{\link{distmap}}.
      \item
      an object of class \code{"nnfun"}
      representing a nearest neighbour function
      (created by the command \code{\link{nnfun}}).
      The fast approximation is \code{\link{nnmap}}.
      \item
      an object of class \code{"densityfun"} representing a kernel
      estimate of intensity (created by
      the command \code{\link[spatstat.explore]{densityfun}}).
      The fast approximation is the Fast Fourier Transform
      algorithm in \code{\link[spatstat.explore]{density.ppp}}.
      \item
      an object of class \code{"Smoothfun"} 
      representing kernel-smoothed values
      (created by the command \code{\link[spatstat.explore]{Smoothfun}}).
      The fast approximation is the Fast Fourier Transform
      algorithm in \code{\link[spatstat.explore]{Smooth.ppp}}.
    }
    \item
    An \code{expression} involving the variables \code{x} and \code{y}
    representing the spatial coordinates, and possibly also involving other
    variables.
    The additional argument \code{W} must be present; it will be converted to
    a binary image mask. The expression \code{X} will be evaluated
    in an environment where \code{x} and \code{y} are
    \bold{vectors} containing the spatial coordinates
    of all the pixels in the image mask.
    Evaluation of the expression \code{X} must yield a vector or factor,
    of the same length as \code{x} and \code{y}, giving the pixel values.
    \item
    a list with entries \code{x, y, z} in the format expected by
    the standard \code{R} functions
    \code{\link{image.default}} and \code{\link{contour.default}}.
    That is, \code{z} is a matrix of pixel values, \code{x} and \code{y}
    are vectors of \eqn{x} and \eqn{y} coordinates respectively,
    and \code{z[i,j]} is the pixel value for the location
    \code{(x[i],y[j])}.
    \item
    a point pattern (object of class \code{"ppp"}).
    See the separate documentation for \code{\link{as.im.ppp}}.
    \item
    A data frame with at least three columns.
    Columns named \code{x}, \code{y} and \code{z}, if present,
    will be assumed to contain the spatial coordinates and the pixel
    values, respectively. Otherwise the \code{x} and \code{y}
    coordinates will be taken from the first two columns of the
    data frame, and any remaining columns will be interpreted as
    pixel values.
  }

  The spatial domain (enclosing rectangle) of the pixel image
  is determined by the argument \code{W}. If \code{W} is absent, 
  the spatial domain is determined by \code{X}.
  When \code{X} is a function, a matrix, or a single numerical value,
  \code{W} is required.
  
  The pixel array dimensions of the final resulting image are determined
  by (in priority order)
  \itemize{
    \item the argument \code{eps}, \code{dimyx} or \code{xy} if present;
    \item the pixel dimensions of the window \code{W}, if it is
    present and if it is a binary mask;
    \item the pixel dimensions of \code{X} if it is an image,
    a binary mask, or a \code{list(x,y,z)};
    \item the default pixel dimensions, 
    controlled by \code{\link{spatstat.options}}.
  }

  Note that if \code{eps}, \code{dimyx} or \code{xy}
  is given, this will override
  the pixel dimensions of \code{X} if it has them.
  Thus, \code{as.im} can be used to change an image's pixel dimensions.

  If the argument \code{na.replace} is given, then all \code{NA} entries
  in the image will be replaced by this value.  The resulting image is
  then defined everwhere on the full rectangular domain, instead of a
  smaller window. Here \code{na.replace} should be a single value,
  of the same type as the other entries in the image.

  If \code{X} is a pixel image that was created by an older version
  of \pkg{spatstat}, the command \code{X <- as.im(X)} will
  repair the internal format of \code{X} so that it conforms to the
  current version of \pkg{spatstat}.

  If \code{X} is a data frame with \code{m} columns,
  then \code{m-2} columns of data are interpreted as pixel values,
  yielding \code{m-2} pixel images. The result of
  \code{as.im.data.frame} is a list of pixel
  images, belonging to the class \code{"imlist"}.
  If \code{m = 3} and \code{drop=TRUE} (the default), then the
  result is a pixel image rather than a list containing this image.

  If \code{X} is a \code{function(x,y)} which returns a matrix of
  values, then \code{as.im(X, W)} will be a list of pixel images.
}
\section{Character-valued images}{
  By default, if the pixel value data are character strings, they will be
  treated as levels of a factor, and the resulting image will be
  factor-valued. To prevent the conversion of character strings to
  factors, use the argument \code{stringsAsFactors=FALSE},
  which is recognised by most of the
  methods for \code{as.im}, or alternatively set
  \code{options(stringsAsFactors=FALSE)}.
}
\section{Handling Character Strings}{
  The argument \code{stringsAsFactors}
  is a logical value (passed to \code{\link[base]{data.frame}})
  specifying how to handle pixel values which
  are character strings. If \code{TRUE}, character values are
  interpreted as factor levels. If \code{FALSE}, they remain
  as character strings. The default values of \code{stringsAsFactors}
  depends on the version of \R.
  \itemize{
    \item In \R versions \code{< 4.1.0}
    the factory-fresh default is \code{stringsAsFactors=FALSE} and
    the default can be changed by setting
    \code{options(stringsAsFactors=FALSE)}.
    \item In \R versions \code{>= 4.1.0} the default
    is \code{stringsAsFactors=FALSE} and there is no option to
    change the default.
  }
}
\value{
  A pixel image (object of class \code{"im"}),
  or a list of pixel images, 
  or \code{NULL} if the conversion failed.
}
\seealso{
  Separate documentation for \code{\link{as.im.ppp}}
}
\examples{
  # window object
  W <- Window(demopat)
  plot(W)
  Z <- as.im(W)
  image(Z)
  # function
  Z <- as.im(function(x,y) {x^2 + y^2}, unit.square())
  image(Z)
  # or as an expression
  Z <- as.im(expression(x^2+y^2), square(1))

  # function with extra arguments
  f <- function(x, y, x0, y0) {
      sqrt((x - x0)^2 + (y-y0)^2)
  }
  Z <- as.im(f, unit.square(), x0=0.5, y0=0.5)
  image(Z)

  # Revisit the Sixties
  Z <- as.im(f, letterR, x0=2.5, y0=2)
  image(Z)
  # usual convention in R
  stuff <- list(x=1:10, y=1:10, z=matrix(1:100, nrow=10))
  Z <- as.im(stuff)
  # convert to finer grid
  Z <- as.im(Z, dimyx=256)

  #' distance functions
  d <- distfun(redwood)
  Zapprox <- as.im(d)
  Zexact <- as.im(d, approx=FALSE)
  plot(solist(approx=Zapprox, exact=Zexact), main="")

  # pixellate the Dirichlet tessellation
  Di <- dirichlet(redwood)
  plot(as.im(Di))
  plot(Di, add=TRUE, border="white")

  # as.im.data.frame is the reverse of as.data.frame.im
  grad <- bei.extra$grad
  slopedata <- as.data.frame(grad)
  slope <- as.im(slopedata)
  unitname(grad) <- unitname(slope) <- unitname(grad) # for compatibility
  all.equal(slope, grad) # TRUE

  ## handling of character values
  as.im("a", W=letterR, na.replace="b")
  as.im("a", W=letterR, na.replace="b", stringsAsFactors=FALSE)
}
\author{
  \spatstatAuthors
}
\keyword{spatial}
\keyword{manip}
