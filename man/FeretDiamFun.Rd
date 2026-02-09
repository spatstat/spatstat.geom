\name{FeretDiamFun}
\alias{FeretDiamFun}
\alias{minFeretDiam}
\alias{maxFeretDiam}
\title{
  Feret Diameter as a Function of Direction
}
\description{
  Create a function which calculates the Feret diameter of a domain
  at any given angle.
  Alternatively, calculate the minimum or maximum value of this function.
}
\usage{
FeretDiamFun(W)
minFeretDiam(W)
maxFeretDiam(W)
}
\arguments{
  \item{W}{
    Window (object of class \code{"owin"}) representing the domain.
  }
}
\details{
  The Feret diameter or caliper diameter of a domain \eqn{W}, at a given angle
  \eqn{\theta}{theta}, is length of the projection of \eqn{W}
  onto a line with orientation \eqn{\theta}{theta}.

  Equivalently, the Feret diameter at angle \eqn{\theta}{theta}
  is the distance between the tangents to \eqn{W}
  that are perpendicular to direction \eqn{\theta}{theta}.

  The Feret diameter at orientation \eqn{\theta}{theta}
  can also be written
  \deqn{
    d_W(\theta) = h_W(\theta) - h_W(\theta + \pi)
  }{
    d[W](theta) = h[W](theta) - h[W](theta + pi)
  }
  where \eqn{h_W(\theta)}{h[W](theta)} is the support function of \code{W}.

  The command \code{FeretDiamFun}
  creates a \code{function} in the \R language, with one argument
  \code{theta}, which can calculate the Feret diameter at any given
  angle or angles \code{theta}.

  The commands \code{minFeretDiam} and \code{maxFeretDiam} calculate
  the minimum and maximum value, respectively, of Feret diameter at any orientation.
}
\value{
  The result of \code{FeretDiamFun} is a
  \code{function} in the \R language with one argument \code{theta}.

  The result of \code{minFeretDiam} or \code{maxFeretDiam} is a single
  numeric value.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{FeretBox}},
  \code{\link{SupportFun}}
}
\examples{
  W <- letterR
  Frame(W) <- grow.rectangle(Frame(W), 0.2)
  plot(W, main="")

  d <- FeretDiamFun(W)
  d(pi/6)

  f <- SupportFun(W)
  tt <- pi/6 + c(0, pi)
  plot(infline(theta=tt, p=f(tt)), col="red", lwd=2)

  ## as a function
  curve(d, to=2*pi)

  minFeretDiam(W)
}
\keyword{spatial}
\keyword{manip}
\keyword{math}
