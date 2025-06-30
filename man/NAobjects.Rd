\name{NAobjects}
\alias{NA_ppp_}
\alias{NA_im_}
\docType{data}
\title{
  Unavailable or Missing Objects
}
\description{
  Objects which represent a missing or unavailable object, belonging to a
  particular class in the \pkg{spatstat} package.
}
\usage{NA_ppp_}
\format{
  \code{NA_ppp_} is an object inheriting the class \code{"ppp"}
  but which is also recognised as \code{NA} by the function \code{\link{is.na}}.

  \code{NA_im_} is an object inheriting the class \code{"im"}
  but which is also recognised as \code{NA} by the function \code{\link{is.na}}.
}
\details{
  Objects which represent a missing or unavailable object, belonging to a
  particular class in the \pkg{spatstat} package.
}
\examples{
  is.na(NA_ppp_)
  is.na(NA_im_)
  X <- list(cells, redwood, cells)
  X[[2]] <- NA_ppp_
  sapply(X, is.na)
}
\keyword{datasets}
