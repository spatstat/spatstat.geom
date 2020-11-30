\name{WindowOnly}
\alias{Window.quad}
\alias{Window.quadratcount}
\alias{Window.tess}
\alias{Window.layered}
\alias{Window.distfun}
\alias{Window.nnfun}
\alias{Window.funxy}
\title{Extract Window of Spatial Object}
\description{
  Given a spatial object (such as a point pattern or pixel image)
  in two dimensions, these functions extract the window
  in which the object is defined.
}
\usage{

 \method{Window}{quad}(X, \dots)

 \method{Window}{quadratcount}(X, \dots)

 \method{Window}{tess}(X, \dots)

 \method{Window}{layered}(X, \dots)

 \method{Window}{distfun}(X, \dots)

 \method{Window}{nnfun}(X, \dots)

 \method{Window}{funxy}(X, \dots)

}
\arguments{
  \item{X}{A spatial object.}
  \item{\dots}{Ignored.}
}
\value{
  An object of class \code{"owin"} (see \code{\link{owin.object}})
  specifying an observation window.
}
\details{
  These are methods for the generic function \code{\link{Window}}
  which extract the spatial window in which the object \code{X}
  is defined. 
}
\seealso{
  \code{\link{Window}},
  \code{\link{Window.ppp}},
  \code{\link{Window.psp}}.
  
  \code{\link{owin.object}}
}
\examples{
   X <- quadratcount(cells, 4)
   Window(X)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{manip}
 
 
