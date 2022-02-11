\name{domain}
\alias{domain}
\alias{domain.ppp}
\alias{domain.psp}
\alias{domain.im}
\alias{domain.ppx}
\alias{domain.pp3}
\alias{domain.quad}
\alias{domain.quadratcount}
\alias{domain.tess}
\alias{domain.layered}
\alias{domain.distfun}
\alias{domain.nnfun}
\alias{domain.funxy}
\title{
  Extract the Domain of any Spatial Object
}
\description{
  Given a spatial object such as a point pattern, in any number of dimensions,
  this function extracts the spatial domain in which the object is defined.
}
\usage{
   domain(X, \dots)

   \method{domain}{ppp}(X, \dots)

   \method{domain}{psp}(X, \dots)

   \method{domain}{im}(X, \dots)

   \method{domain}{ppx}(X, \dots)

   \method{domain}{pp3}(X, \dots)

 \method{domain}{quad}(X, \dots)

 \method{domain}{quadratcount}(X, \dots)

 \method{domain}{tess}(X, \dots)

 \method{domain}{layered}(X, \dots)

 \method{domain}{distfun}(X, \dots)

 \method{domain}{nnfun}(X, \dots)

 \method{domain}{funxy}(X, \dots)




}
\arguments{
  \item{X}{
    A spatial object such as a point pattern (in any number
    of dimensions), line segment pattern or pixel image.
  }
  \item{\dots}{
    Extra arguments. They are ignored by all the methods listed here.
  }
}
\details{
  The function \code{\link[spatstat.geom]{domain}} is generic.

  For a spatial object \code{X} in any number of dimensions, 
  \code{domain(X)} extracts the spatial domain in which \code{X} is
  defined.

  For a two-dimensional object \code{X}, typically \code{domain(X)}
  is the same as \code{Window(X)}.

  Exceptions occur for methods related to linear networks.

}
\value{
  A spatial object representing the domain of \code{X}.
  Typically a window (object of class \code{"owin"}),
  a three-dimensional box (\code{"box3"}), a multidimensional
  box (\code{"boxx"}) or a linear network (\code{"linnet"}).
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link[spatstat.core]{domain.ppm}},
  \code{\link[spatstat.random]{domain.rmhmodel}},
  \code{\link[spatstat.linnet]{domain.lpp}}.
  \code{\link{Window}},
  \code{\link{Frame}}.
}
\examples{
  domain(redwood)
  domain(bei.extra$elev)
  domain(osteo$pts[[1]])
}
\keyword{spatial}
\keyword{manip}
