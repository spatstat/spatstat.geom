\name{default.image.colours}
\alias{default.image.colours}
\alias{reset.default.image.colours}
\title{
  Default Colours for Images in Spatstat 
}
\description{
  Extract or change the default colours for images in \pkg{spatstat}.
}
\usage{
default.image.colours()
reset.default.image.colours(col = NULL)
}
\arguments{
  \item{col}{
    A vector of colour values.
  }
}
\details{
  These functions extract and change the current default colours
  used for plotting colour images in the \pkg{spatstat} family of
  packages, in particular by the functions
  \code{\link[spatstat.geom]{plot.im}}
  and \code{\link[spatstat.linnet]{plot.linim}}.

  The default colour values are a vector of
  character strings which can be interpreted as colours.
  In any particular instance of
  \code{\link[spatstat.geom]{plot.im}}
  or \code{\link[spatstat.linnet]{plot.linim}},
  the default colours are interpolated to obtain a vector of colour
  values of the required length (usually 256, controlled by the
  argument \code{ncolours} to the plot command).
  
  \code{default.image.colours()} returns the current default colours.
  \code{reset.default.image.colours(col)} sets the default colours to
  be the vector \code{col}.
  \code{reset.default.image.colours()} or
  \code{reset.default.image.colours(NULL)} resets the
  factory default, which is row 29 of the Kovesi uniform
  perceptual contrast table described in \code{\link[spatstat.data]{Kovesi}}.
}
\value{
  A character vector of values which can be interpreted as colours.
}
\author{
  \adrian.
}
\seealso{
  \code{\link[spatstat.geom]{plot.im}}
}
\examples{
 a <- default.image.colours()
 length(a)
}
\keyword{color}
\keyword{hplot}
