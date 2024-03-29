\name{anylist}
\alias{anylist}
\alias{as.anylist}
\title{
  List of Objects
}
\description{
  Make a list of objects of any type.
}
\usage{
anylist(\dots)
as.anylist(x)
}
\arguments{
  \item{\dots}{
    Any number of arguments of any type.
  }
  \item{x}{
    A list.
  }
}
\details{
  An object of class \code{"anylist"} is
  a list of objects that the user intends to treat in a similar fashion.

  For example it may be desired to plot each of the objects
  side-by-side: this can be done using the function
  \code{\link{plot.anylist}}.
  
  The objects can belong to any class;
  they may or may not all belong to the same class.

  In the \pkg{spatstat} package, various functions produce
  an object of class \code{"anylist"}.
}
\value{
  A list, belonging to the class \code{"anylist"},
  containing the original objects.
}
\author{\adrian
  
  
  \rolf
  
  and \ege
  
}
\seealso{
  \code{\link{solist}},
  \code{\link{as.solist}},
  \code{\link{anylapply}}.
}
\examples{
  if(require(spatstat.explore)) {
    anylist(cells, intensity(cells), Kest(cells))
  } else {
    anylist(cells, intensity(cells))
  }
  anylist()
}
\keyword{list}
\keyword{manip}
