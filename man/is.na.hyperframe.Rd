\name{is.na.hyperframe}
\alias{is.na.hyperframe}
\title{
  Identify Missing Entries in a Hyperframe
}
\description{
  Given a hyperframe, this function returns a logical matrix
  specifying which entries of the hyperframe are missing.
}
\usage{
\method{is.na}{hyperframe}(x)
}
\arguments{
  \item{x}{
    A hyperframe (object of class \code{"hyperframe"}).
  }
}
\details{
  This function returns a logical matrix, with the same dimensions
  as the hyperframe \code{x}, with the value \code{TRUE}
  in each position where the corresponding entry in the hyperframe
  is missing.

  An entry \code{x[i,j]} is deemed to be missing if it is
  either \code{NA} (representing a missing value in an atomic vector)
  or is an object of class \code{"NAobject"} (representing a missing
  object of a particular class in the \pkg{spatstat} package family). 
}
\value{
  A logical matrix with the same dimensions as \code{x}.
}
\author{
  \adrian and \ege.
}
\seealso{
  \code{\link{NAobject}}
}
\examples{
  A <- c(1, NA, 3)
  B <- solist(cells, cells, NA_ppp_)
  h <- hyperframe(A=A, B=B)
  h
  is.na(h)
}
\keyword{spatial}
\keyword{manip}
