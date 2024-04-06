\name{default.symbolmap}
\alias{default.symbolmap}
\title{
  Default Symbol Map for Plotting a Spatial Pattern
}
\description{
  Determines the symbol map for plotting a spatial pattern,
  when one is not supplied by the user.
}
\usage{
default.symbolmap(x, \dots)
}
\arguments{
  \item{x}{
    A spatial object in the \pkg{spatstat} package,
    such as a point pattern (class \code{"ppp"}).
  }
  \item{\dots}{
    Additional arguments passed to methods.
  }
}
\details{
  In the \pkg{spatstat} package, an object of class \code{"symbolmap"}
  defines a mapping between data and graphical symbols.

  If a plot command \code{plot(x, \dots)} has been issued,
  and if the arguments were not sufficient to determine the
  symbol map that should be used, then \code{default.symbolmap(x, \dots)} 
  will be executed to determine the default symbol map.

  The function \code{default.symbolmap} is generic, with a method
  for point patterns (class \code{"ppp"}) and possibly for other
  classes.
}
\value{
  A symbol map (object of class \code{"symbolmap"}).
}
\author{
  \adrian
}
\seealso{
  \code{\link{default.symbolmap.ppp}}
}
\keyword{spatial}
\keyword{hplot}
