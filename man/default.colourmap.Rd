\name{default.colourmap}
\alias{default.colourmap}
\title{
  Default Colour Map for Plotting a Spatial Pattern
}
\description{
  Determines a colour map for plotting a spatial pattern
  (or other data) when one is not supplied by the user.
}
\usage{
default.colourmap(x, \dots,
                  col=spatstat.options("image.colfun"),
                  scramble.cols=FALSE,
                  monochrome=spatstat.options("monochrome"))
}
\arguments{
  \item{x}{
    A vector of atomic values, or a factor.
  }
  \item{col}{
    Optional. A specification of colours to be used. See Details.
  }
  \item{scramble.cols}{
    Logical value. If \code{TRUE}, the sequence of colour values
    will be randomly permuted. 
  }
  \item{monochrome}{
    Logical value. If \code{TRUE}, the colours will be converted
    to greyscale colours.
  }
  \item{\dots}{
    Additional arguments passed to methods.
  }
}
\details{
  In the \pkg{spatstat} package, an object of class \code{"colourmap"}
  defines a mapping between data and colours.

  The function \code{default.colourmap} provides a suitable default
  colour map for the values in \code{x}.

  If \code{x} is a factor, \code{default.colourmap(x)} is a mapping
  from the factor levels to colours.
  If \code{x} is a logical vector, \code{default.colourmap(x)} is a mapping
  from the values \code{TRUE} and \code{FALSE} to colours.
  If \code{x} is a numeric vector, \code{default.colourmap(x)} is a
  mapping from numbers in the interval between the minimum and maximum values of
  \code{x} to colours.

  The argument \code{col} may provide colour information
  in any of the following formats:
  \itemize{
    \item A colour map (object of class \code{"colourmap"})
    \item A palette function (a \code{function(n)} or \code{function(n, \dots)} 
    which returns a vector of colour values,
    such as \code{\link[grDevices]{rainbow}}
    \item A vector of integers between 1 and 8 indexing the standard
    colour palette
    \item A vector of character strings giving common names of colours
    (e.g. \code{"red"})
    \item A vector of character strings giving hexadecimal codes for
    colours (e.g. produced by \code{\link[grDevices]{rgb}}
    or \code{\link[grDevices]{rainbow}})
    \item \code{NULL}.
  }
}
\section{Scrambling the colours}{
  The argument \code{scramble.cols=TRUE}
  will cause the sequence of colours in the colour map to be
  randomly permuted. This is a useful trick when it is desired that
  adjacent colours in the sequence should be easily distinguishable.
  The return value will be the randomised colour map.

  Randomisation implies that the colour map will be different
  each time the command is executed. For reproducible results,
  set the random generator seed using \code{\link{set.seed}}.
}
\value{
  A colour map (object of class \code{"colourmap"}).
}
\author{
  \adrian
}
\seealso{
  \code{\link{colourmap}}
}
\examples{
  default.colourmap(letters[1:4])
  default.colourmap(factor(letters[1:3], levels=letters[1:4]))
  default.colourmap(TRUE)
  default.colourmap(NULL)
}
\keyword{spatial}
\keyword{hplot}
