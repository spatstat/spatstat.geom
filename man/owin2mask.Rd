\name{owin2mask}
\alias{owin2mask}
\title{
  Convert Window to Binary Mask under Constraints 
}
\description{
  Converts a window (object of class \code{"owin"})
  to a binary pixel mask, with options to require
  that the inside, outside, and/or boundary of the window
  should be completely covered.
}
\usage{
owin2mask(W,
         op = c("sample", "notsample",
                "cover", "inside",
                "uncover", "outside", "boundary"),
         \dots)
}
\arguments{
  \item{W}{
    A window (object of class \code{"owin"}).
  }
  \item{op}{
    Character string (partially matched) specifying how \code{W}
    should be converted to a binary pixel mask.
  }
  \item{\dots}{
    Optional arguments passed to \code{\link{as.mask}} to determine
    the pixel resolution.
  }
}
\details{
  This function is similar to, but more flexible than, \code{\link{as.mask}}.
  It converts the interior, exterior, or boundary of the window \code{W}
  to a binary pixel mask.

  \itemize{
    \item
    If \code{op="sample"} (the default), the mask consists of all pixels
    whose \bold{centres} lie inside the window \code{W}.
    This is the same as using \code{\link{as.mask}}.
    \item 
    If \code{op="notsample"}, the mask consists of all pixels
    whose \emph{centres lie outside} the window \code{W}.
    This is the same as using \code{\link{as.mask}} followed by
    \code{\link{complement.owin}}.
    \item 
    If \code{op="cover"}, the mask consists of all pixels
    which overlap the window \code{W}, wholly or partially. 
    \item 
    If \code{op="inside"}, the mask consists of all pixels which
    lie entirely inside the window \code{W}.
    \item 
    If \code{op="uncover"}, the mask consists of all pixels which
    overlap the outside of the window \code{W}, wholly or partially.
    \item 
    If \code{op="outside"}, the mask consists of all pixels which
    lie entirely outside the window \code{W}.
    \item 
    If \code{op="boundary"}, the mask consists of all pixels which
    overlap the boundary of the window \code{W}.
  }
}
\value{
  A mask (object of class \code{"owin"} of type \code{"mask"}
  representing a binary pixel mask).
}
\author{
  \adrian.
}
\seealso{
  \code{\link{as.mask}}
}
\examples{
  W <- Window(chorley)
  opa <- par(mfrow=c(2,4))
  plot(as.mask(W, dimyx=16), col="grey", main="as.mask")
  plot(W, add=TRUE, lwd=3, border="red")  
  plot(owin2mask(W, "i", dimyx=16), col="grey", main="inside")
  plot(W, add=TRUE, lwd=3, border="red")
  plot(owin2mask(W, "c", dimyx=16), col="grey", main="cover")
  plot(W, add=TRUE, lwd=3, border="red")
  plot(owin2mask(W, "b", dimyx=16), col="grey", main="boundary")
  plot(W, add=TRUE, lwd=3, border="red")
  plot(as.mask(complement.owin(W), dimyx=16), col="grey",
       main="as.mask(complement.owin)")
  plot(W, add=TRUE, lwd=3, border="red")  
  plot(owin2mask(W, "o", dimyx=16), col="grey", main="outside")
  plot(W, add=TRUE, lwd=3, border="red")
  plot(owin2mask(W, "u", dimyx=16), col="grey", main="uncover")
  plot(W, add=TRUE, lwd=3, border="red")
  plot(owin2mask(W, "b", dimyx=16), col="grey", main="boundary")
  plot(W, add=TRUE, lwd=3, border="red")
  par(opa)
}
\keyword{spatial}
\keyword{manip}
