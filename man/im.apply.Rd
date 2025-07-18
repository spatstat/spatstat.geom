\name{im.apply}
\alias{im.apply}
\title{
  Apply Function Pixelwise to List of Images 
}
\description{
  Returns a pixel image obtained by applying a given function
  to corresponding pixels in several pixel images.
}
\usage{
im.apply(X, FUN, \dots, fun.handles.na=FALSE, check=TRUE, verbose=TRUE)
}
\arguments{
  \item{X}{
    A list of pixel images (objects of class \code{"im"}).
  }
  \item{FUN}{
    A function that can be applied to vectors,
    or a character string giving the name of such a function.
  }
  \item{\dots}{
    Additional arguments to \code{FUN}.
  }
  \item{fun.handles.na}{
    Logical value specifying what to do when the data include
    \code{NA} values. See Details.
  }
  \item{check}{
    Logical value specifying whether to check that the images in
    \code{X} are compatible (for example that they have the same
    grid of pixel locations) and to convert them to compatible images
    if necessary.
  }
  \item{verbose}{
    Logical value specifying whether to print informative messages.
  }
}
\details{
  The argument \code{X} should be a list of pixel images
  (objects of class \code{"im"}). If the images do not have
  identical pixel grids, they will be converted to a common
  grid using \code{\link{harmonise.im}}.

  At each pixel location, the values of the images in \code{X}
  at that pixel will be extracted as a vector; the function
  \code{FUN} will be applied to this vector; and the return value of
  \code{FUN} will become the pixel value of the resulting image.
  For example \code{im.apply(X, mean)} will return a pixel image
  in which the value of each pixel is the average of the corresponding
  pixel values in the images in \code{X}.

  If the result of \code{FUN} is a vector, then the result of
  \code{im.apply} will be a list of images. For example \code{im.apply(X, range)} will
  return a list of two images containing the pixelwise minimum
  and pixelwise maximum, respectively, of the input images in \code{X}.

  The argument \code{fun.handles.na} specifies what to do
  when some of the pixel values are \code{NA}.
  \itemize{
    \item If \code{fun.handles.na=FALSE} (the default),
    the function \code{FUN} is never applied to data that include
    \code{NA} values; the result is defined to be \code{NA} whenever
    the data contain \code{NA}.
    \item
    If \code{fun.handles.na=TRUE}, the function \code{FUN} will be applied to
    all pixel data, including those which contain \code{NA} values.
  }
}
\section{Computation details}{
  The computation is performed efficiently
  using matrix operations where possible,
  by constructing a large matrix containing all the pixel data for all images.
  However, if this matrix would exceed the maximum permitted size of matrices
  (as given in \code{\link{spatstat.options}("maxmatrix")}),
  then the calculation is performed on smaller sub-images,
  and this is reported in a message if \code{verbose=TRUE}.
}
\value{
  A pixel image (object of class \code{"im"}) or a list of pixel images.
}
\seealso{
  \code{\link{eval.im}} for algebraic operations with images.
}
\examples{
  # list of two pixel images
  Y <- solapply(bei.extra, scaletointerval)
  plot(Y)
  im.apply(Y, max)
  im.apply(Y, sum)
  im.apply(Y, range)

  ## Example with incompatible patterns of NA values
  B <- owin(c(438, 666), c(80, 310))
  Y[[1]][B] <- NA
  opa <- par(mfrow=c(2,2))
  plot(Y[[1]])
  plot(Y[[2]])
  #' Default action: NA -> NA
  plot(im.apply(Y, mean))
  #' Use NA handling in mean.default
  plot(im.apply(Y, mean, na.rm=TRUE, fun.handles.na=TRUE))
  par(opa)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{manip}
\keyword{programming}
