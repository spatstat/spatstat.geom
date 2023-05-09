\name{is.boxx}
\Rdversion{1.1}
\alias{is.boxx}
\title{
  Recognise a Multi-Dimensional Box
}
\description{
  Checks whether its argument is a multidimensional box
  (object of class \code{"boxx"}).

}
\usage{
  is.boxx(x)
}
\arguments{
  \item{x}{
    Any object.
  }
}
\details{
  This function tests whether the object \code{x}
  is a multidimensional box of class \code{"boxx"}.

  The result is determined to be \code{TRUE} if \code{x}
  inherits from \code{"boxx"}, i.e. if \code{x}
  has \code{"boxx"} amongst its classes.
}
\value{
  A logical value.
}
\author{
  \spatstatAuthors
}
\seealso{
  \code{\link{methods.boxx}},
  \code{\link{boxx}}.
}
\examples{
   B <- boxx(c(0,10),c(0,10),c(0,5),c(0,1), unitname="km")
   is.boxx(B)
   is.boxx(42)
}
\keyword{spatial}
\keyword{methods}
