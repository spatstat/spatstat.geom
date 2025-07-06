\name{is.NAobject}
\alias{is.NAobject}
\title{
  Recognise NA Objects
}
\description{
  Recognises whether an object is
  an NA object (representing a missing or unavailable object
  in the \pkg{spatstat} package family).
}
\usage{
  is.NAobject(x)
}
\arguments{
  \item{x}{
    An object.
  }
}
\details{
  This function recognises whether an object is
  a missing or unavailable object belonging to one of the
  classes in the \pkg{spatstat} package family.

  In \pkg{spatstat}, a missing or unavailable object of class
  \code{"foo"} is represented by an object that inherits
  the classes \code{"foo"} and \code{"NAobject"}. 

  For any object \code{x}, the command
  \code{is.NAobject(x)} will return \code{TRUE} if
  \code{x} is a missing or unavailable object, and \code{FALSE} otherwise.
}
\value{
  A single logical value.
}
\author{
  \adrian and \ege.
}
\seealso{
  \code{\link{NAobject}}
}
\examples{
   A <- NAobject("ppp")
   is.NAobject(A)

   X <- solist(cells, NA, cells)
   is.NAobject(X[[2]])
}
\keyword{methods}
