\name{fillholes.owin}
\Rdversion{1.1}
\alias{fillholes.owin}
\title{
  Remove Small Holes in a Window
}
\description{
  Given a window, this function removes any small holes in the window.
}
\usage{
fillholes.owin(W, amin)
}
\arguments{
  \item{W}{
    Window (object of class \code{"owin"}).
  }
  \item{amin}{
    Numeric value. The smallest permissible area of a hole.
  }
}
\details{
  This function simplifies a window \code{W}
  by removing any holes with area less than \code{amin}.

  The argument \code{W} must be a window (object of class
  \code{"owin"}). 
}
\value{
  Another window (object of class \code{"owin"})
  of the same type as \code{W}.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{simplify.owin}},
  \code{\link{owin}}
}
\examples{
  ## make object with small hole
  A <- meningitis$kreise
  ta <- tile.areas(A)
  smalltile <- tiles(A)[[which.min(ta)]]
  W <- setminus.owin(Window(A), smalltile)
  amin <- sort(ta)[3]

  opa <- par(mfrow=c(1,2))
  plot(W)
  plot(fillholes.owin(W, amin))

  M <- as.mask(W)
  plot(M)
  plot(fillholes.owin(M, amin))
}
\keyword{spatial}
\keyword{math}
