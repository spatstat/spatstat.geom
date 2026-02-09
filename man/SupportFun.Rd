\name{SupportFun}
\alias{SupportFun}
\title{
  Support Function
}
\description{
  Create the support function of a two-dimensional domain.
}
\usage{
SupportFun(W, origin = c(0, 0))
}
\arguments{
  \item{W}{
    Window (object of class \code{"owin"}) representing the domain.
  }
  \item{origin}{
    The location of the origin for the support function.
    Either a numeric vector of length 2 giving the location,
    or a point pattern containing only one point,
    or a list with two entries named \code{x} and \code{y},
    or one of the character strings
    \code{"centroid"}, \code{"midpoint"},
    \code{"left"}, \code{"right"}, \code{"top"}, \code{"bottom"},
    \code{"topleft"}, \code{"bottomleft"}, \code{"topright"} or
    \code{"bottomright"} (partially matched).
  }
}
\details{
  The support function \eqn{h_W(\theta)}{h[W](theta)} of a set \eqn{W}
  is defined by
  \deqn{
    h_W(\theta) = \sup\{ x \cos\theta + y\sin\theta : (x,y) \in W \}
  }{
    h[W](theta) = sup( x * cos(theta) + y * sin(theta): (x,y) in W)
  }
  defined for all angles \eqn{theta} in \eqn{(0, 2\pi)}{(0, 2*pi)}.

  For any direction \code{theta} (in radians),
  \eqn{h_W(\theta)}{h[W](theta)} is the signed distance from the origin
  to a line which is tangent to \code{W} and perpendicular to \code{theta}.

  The result of \code{SupportFun(W)} is a \code{function} in the \R
  language with one argument \code{theta}. This function can be applied
  to any vector of values of \eqn{\theta}{theta} to compute the values
  of the support function.
}
\value{
  A \code{function} in the \R language with one argument \code{theta}.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{FeretDiamFun}},
  \code{\link{polarbody}}
}
\examples{
  W <- letterR
  Frame(W) <- grow.rectangle(Frame(W), 0.2)
  plot(W, main="")

  f <- SupportFun(W)
  p <- f(pi/6)
  plot(infline(p=p, theta=pi/6), col="red", lwd=2)

  ## as a function
  curve(f, to=2*pi)
}
\keyword{spatial}
\keyword{manip}
\keyword{math}
