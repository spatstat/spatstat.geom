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
                "uncover", "outside",
                "boundary",
                "majority", "minority"),
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
    Optional arguments passed to \code{\link[spatstat.geom]{as.mask}} to determine
    the pixel resolution.
  }
}
\details{
  This function is similar to, but more flexible than, \code{\link[spatstat.geom]{as.mask}}.
  It converts the interior, exterior, or boundary of the window \code{W}
  to a binary pixel mask.

  \itemize{
    \item
    If \code{op="sample"} (the default), the mask consists of all pixels
    whose \bold{centres} lie inside the window \code{W}.
    This is the same as using \code{\link[spatstat.geom]{as.mask}}.
    \item 
    If \code{op="notsample"}, the mask consists of all pixels
    whose \emph{centres lie outside} the window \code{W}.
    This is the same as using \code{\link[spatstat.geom]{as.mask}} followed by
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
    \item 
    If \code{op="majority"}, the mask consists of all pixels in which
    at least half of the pixel area is covered by the window \code{W}.
    \item 
    If \code{op="minority"}, the mask consists of all pixels in which
    less than half of the pixel area is covered by the window \code{W}.
  }
  These operations are complementary pairs as follows:
  \tabular{lll}{
    \code{"notsample"} \tab is the complement of \tab \code{"sample"} \cr
    \code{"uncover"} \tab is the complement of \tab \code{"inside"} \cr
    \code{"outside"} \tab is the complement of \tab \code{"cover"} \cr
    \code{"minority"} \tab is the complement of \tab \code{"majority"}
  }
  They also satisfy the following set relations:
  \tabular{lll}{
    \code{"inside"} \tab is a subset of \tab \code{"cover"} \cr
    \code{"outside"} \tab is a subset of \tab \code{"uncover"} \cr
    \code{"boundary"} \tab is a subset of \tab \code{"cover"} \cr
    \code{"boundary"} \tab is a subset of \tab \code{"uncover"}
  }
  The results of \code{"inside"}, \code{"boundary"} and \code{"outside"}
  are disjoint and their union is the entire frame.

  Theoretically \code{"sample"} should be a subset of \code{"cover"},
  \code{"notsample"} should be a subset of \code{"uncover"},
  \code{"inside"} should be a subset of \code{"majority"}
  and \code{"outside"} should be a subset of \code{"minority"},
  but these could be false due to numerical error in computational geometry.
}
\value{
  A mask (object of class \code{"owin"} of type \code{"mask"}
  representing a binary pixel mask).
}
\author{
  \adrian.
}
\seealso{
  \code{\link[spatstat.geom]{as.mask}}
}
\examples{
  W <- Window(chorley)
  opa <- par(mfrow=c(2,5))
  plot(as.mask(W, dimyx=10), col="grey", main="sample")
  plot(W, add=TRUE, lwd=3, border="red")  
  plot(owin2mask(W, "ma", dimyx=10), col="grey", main="majority")
  plot(W, add=TRUE, lwd=3, border="red")
  plot(owin2mask(W, "i", dimyx=10), col="grey", main="inside")
  plot(W, add=TRUE, lwd=3, border="red")
  plot(owin2mask(W, "c", dimyx=10), col="grey", main="cover")
  plot(W, add=TRUE, lwd=3, border="red")
  plot(owin2mask(W, "b", dimyx=10), col="grey", main="boundary")
  plot(W, add=TRUE, lwd=3, border="red")
  plot(as.mask(complement.owin(W), dimyx=10), col="grey",
       main="notsample")
  plot(W, add=TRUE, lwd=3, border="red")  
  plot(owin2mask(W, "mi", dimyx=10), col="grey", main="minority")
  plot(W, add=TRUE, lwd=3, border="red")
  plot(owin2mask(W, "o", dimyx=10), col="grey", main="outside")
  plot(W, add=TRUE, lwd=3, border="red")
  plot(owin2mask(W, "u", dimyx=10), col="grey", main="uncover")
  plot(W, add=TRUE, lwd=3, border="red")
  plot(owin2mask(W, "b", dimyx=10), col="grey", main="boundary")
  plot(W, add=TRUE, lwd=3, border="red")
  par(opa)
}
\keyword{spatial}
\keyword{manip}
