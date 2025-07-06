\name{is.na.solist}
\alias{is.na.solist}
\title{
  Identify Missing Entries in a List of Spatial Objects
}
\description{
  Given a list of spatial objects (of class \code{"solist"}),
  this function returns a logical vector
  specifying which entries of the list are missing.
}
\usage{
\method{is.na}{solist}(x)
}
\arguments{
  \item{x}{
    A list of spatial objects (an object of class \code{"solist"}).
  }
}
\details{
  This function returns a logical vector with the same length
  as the list \code{x}, with the value \code{TRUE}
  in each position where the corresponding entry in the list
  is missing.

  An entry \code{x[[i]]} is deemed to be missing if it is
  an object of class \code{"NAobject"} (representing a missing
  object of a particular class in the \pkg{spatstat} package family). 
}
\value{
  A logical vector with the same length as \code{x}.
}
\author{
  \adrian and \ege.
}
\seealso{
  \code{\link{NAobject}}
}
\examples{
  B <- solist(cells, NA, cells)
  B[[2]]
  is.na(B)
  V <- solist(cells, Window(cells), NAobject("im"))
  is.na(V)
}
\keyword{spatial}
\keyword{manip}
