\name{spatstat.geom-deprecated}
\alias{spatstat.geom-deprecated} %DoNotExport
\alias{as.psp.owin}
\alias{circumradius}
\alias{circumradius.owin}
\alias{circumradius.ppp}
\title{Deprecated spatstat.geom functions}
\description{
  Deprecated spatstat.geom functions.
}
\usage{
\method{as.psp}{owin}(x, \dots, window=NULL, 
       check=spatstat.options("checksegments"), fatal=TRUE)
circumradius(x, \dots)
\method{circumradius}{owin}(x, \dots)
\method{circumradius}{ppp}(x, \dots)
}

\details{
  These functions are deprecated, and will eventually be deleted from
  the \pkg{spatstat.geom} package.
  
  \code{as.psp.owin} has been replaced by \code{\link{edges}}.

  \code{circumradius} is replaced by the more appropriately named
  \code{boundingradius}.

  
}
\value{
\code{as.psp.owin} returns a window (object of class \code{"owin"}).

\code{circumradius} returns a numeric value.
}
\keyword{internal}
