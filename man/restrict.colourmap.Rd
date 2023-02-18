\name{restrict.colourmap}
\alias{restrict.colourmap}
\title{
  Restrict a Colour Map to a Subset of Values
}
\description{
  Given a colour map defined on a range of numerical values or a set of
  discrete inputs, the command restricts the range of values to a narrower range,
  or restricts the set of inputs to a subset, and returns the associated
  colour map.
}
\usage{
restrict.colourmap(x, ..., range = NULL, breaks = NULL, inputs = NULL)
}
\arguments{
  \item{x}{
    Colour map (object of class \code{"colourmap"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{range}{
    New, restricted range of numerical values to which the
    colour map will apply. A numeric vector of length 2 giving the
    minimum and maximum values of the input.
    Incompatible with \code{breaks} and \code{inputs}.
  }
  \item{breaks}{
    Vector of breakpoints for the new colour map.
    A numeric vector with increasing entries.
    Incompatible with \code{range} and \code{inputs}.
  }
  \item{inputs}{
    Values accepted as inputs for the new colour map.
    A factor or vector.
    Incompatible with \code{breaks} and \code{range}.
  }
}
\details{
  This command produces a new colour map \code{y} which is consistent with the
  original colour map \code{x}, except that \code{y} is defined on a
  narrower interval of numeric values, or a smaller set of discrete
  input values, than \code{x}.
}
\value{
  Colour map (object of class \code{"colourmap"}).
}
\author{
  \adrian.
}
\seealso{
  \code{\link{colourmap}}
}
\examples{
  plot(a <- colourmap(topo.colors(128), range=c(-1,1)))
  plot(b <- restrict.colourmap(a, range=c(0,1)))
}
\keyword{spatial}
\keyword{color}
