\name{as.colourmap}
\alias{as.colourmap}
\alias{as.colourmap.colourmap}
\alias{as.colourmap.symbolmap}
\title{
  Convert to Colour Map
}
\description{
  Convert some other kind of data to a colour map.
}
\usage{
  as.colourmap(x, \dots)

  \method{as.colourmap}{colourmap}(x, \dots)

  \method{as.colourmap}{symbolmap}(x, \dots, warn=TRUE)
}
\arguments{
  \item{x}{
    Data to be converted to a colour map.
    An object of class \code{"symbolmap"}, \code{"colourmap"}
    or some other kind of suitable data.
  }
  \item{\dots}{
    Other arguments passed to methods.
  }
  \item{warn}{
    Logical value specifying whether to issue a warning
    if \code{x} does not contain any colour map information.
  }
}
\details{
  If \code{x} contains colour map information, it will be extracted
  and returned as a colour map object. Otherwise, \code{NULL} will be
  returned (and a warning will be issued if \code{warn=TRUE}, the default).
}
\value{
  A colour map (object of class \code{"colourmap"}) or \code{NULL}.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{colourmap}}
}
\examples{
  m <- pHcolourmap(c(3,8))
  g <- symbolmap(pch=21, bg=m, size=function(x){ 1.1 * x }, range=c(3,8))
  opa <- par(mfrow=c(1,2))
  plot(g, vertical=TRUE)
  plot(as.colourmap(g), vertical=TRUE)
  par(opa)
}
\keyword{spatial}
\keyword{color}

