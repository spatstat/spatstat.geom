\name{rev.colourmap}
\alias{rev.colourmap}
\title{
  Reverse the Colours in a Colour Map
}
\description{
  Reverse the sequence of colour values in a colour map.
}
\usage{
  \method{rev}{colourmap}(x)
}
\arguments{
  \item{x}{
    A colour map (object of class \code{"colourmap"}).
  }
}
\details{
  This is a method for the generic \code{\link[base]{rev}}
  for the class of colour maps.

  The sequence of colour values in the colour map will be reversed,
  without changing any other details. 
}
\value{
  A colour map (object of class \code{"colourmap"}).
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{colourmap}},
  \code{\link{colouroutputs}}
}
\examples{
  co <- colourmap(rainbow(100), range=c(-1,1))
  opa <- par(mfrow=c(1,2))
  plot(co, vertical=TRUE)
  plot(rev(co), vertical=TRUE)
  par(opa)
}
\keyword{spatial}
\keyword{color}

