\name{is.NAobject}
\alias{is.NAobject}
\alias{is.na.NAobject}
\alias{is.na.ppp}
\alias{is.na.im}
\alias{is.na.owin}
\title{
  Recognise NA Objects
}
\description{
  Functions which recognise whether an object is
  an NA object (representing a missing or unavailable object
  in the \pkg{spatstat} package family).
}
\usage{
  is.NAobject(x)

  \method{is.na}{ppp}(x)
  \method{is.na}{im}(x)
  \method{is.na}{owin}(x)

  \method{is.na}{NAobject}(x)
}
\arguments{
  \item{x}{
    An object.
  }
}
\details{
  These functions recognise whether an object is
  a missing or unavailable object belonging to one of the
  classes in the \pkg{spatstat} package family.

  In \pkg{spatstat}, a missing or unavailable object of class
  \code{"foo"} is represented by an object that inherits
  the classes \code{"foo"} and \code{"NAobject"}. 
  For example, \code{\link{NA_ppp_}} is an object of class \code{"ppp"}
  which also inherits the class \code{"NAobject"} and represents
  a point pattern object which is missing or unavailable.

  For an object \code{x} of any class, the command
  \code{is.NAobject(x)} will return \code{TRUE} if
  \code{x} is a missing or unavailable object, and \code{FALSE} otherwise.

  If an object \code{x} is known to belong to class \code{"ppp"} (for
  example), then \code{is.na(x)} is dispatched to
  the method \code{is.na.ppp} and will return \code{TRUE} if
  \code{x} is missing or unavailable, and \code{FALSE} otherwise.
}
\value{
  A single logical value.
}
\author{
  \adrian and \ege.
}
\seealso{
  \code{\link{NA_ppp_}}
}
\examples{
   X <- solist(cells, NA_ppp_, cells)
   sapply(X, is.na)

   is.NAobject(NA_im_)
}
\keyword{methods}
