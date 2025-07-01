\name{NAobjects}
\alias{NA_ppp_}
\alias{NA_im_}
\alias{NA_owin_}
\docType{data}
\title{
  Unavailable or Missing Objects
}
\description{
  Objects which represent a missing or unavailable object, belonging to a
  particular class in the \pkg{spatstat} package.
}
\usage{
NA_ppp_
NA_im_
NA_owin_
}
\format{
  \code{NA_ppp_},
  \code{NA_im_}
  and \code{NA_owin_}
  are objects inheriting the classes \code{"ppp"}, \code{"im"} and
  \code{"owin"} respectively. They also inherit the class
  \code{"NAobject"}.
}
\details{
  Each of these objects represents a missing or unavailable object
  of a particular class in the \pkg{spatstat} package:
  \itemize{
    \item \code{NA_ppp_} belongs to the class \code{"ppp"}
    \item \code{NA_im_} belongs to the class \code{"im"}
    \item \code{NA_owin_} belongs to the class \code{"owin"}
  }
  Additionally each object belongs to the class \code{"NAobject"}.
}
\section{Recognising an NA object}{
  For any object of class \code{"ppp"}, \code{"im"} or \code{"owin"},
  the generic \code{\link[base]{is.na}} will return \code{TRUE} if the object
  is an NA object, and \code{FALSE} otherwise. This is implemented by
  the methods \code{\link{is.na.ppp}}, \code{\link{is.na.im}} and
  \code{\link{is.na.owin}}. 

  If you are not sure of the class of an object, use
  \code{\link{is.NAobject}} or \code{\link[base]{inherits}(, what="NAobject")}
  to determine whether it is an NA object.
}
\section{Warning}{
  These objects are experimental and are not supported by most of the
  code in the \pkg{spatstat} package family.
}
\seealso{
  \code{\link{is.NAobject}},
  \code{\link{is.na.ppp}},
  \code{\link{is.na.im}}
  \code{\link{is.na.owin}}
}
\examples{
  is.na(NA_ppp_)
  is.na(NA_im_)
  X <- list(cells, redwood, cells)
  X[[2]] <- NA_ppp_
  sapply(X, is.na)
}
\keyword{datasets}
