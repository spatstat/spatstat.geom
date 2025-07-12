\name{methods.layered} 
\Rdversion{1.1}
\alias{methods.layered} %DoNotExport
\alias{c.layered}
\alias{shift.layered}
\alias{reflect.layered}
\alias{flipxy.layered}
\alias{rotate.layered}
\alias{affine.layered}
\alias{rescale.layered}
\alias{scalardilate.layered}
\title{
  Methods for Layered Objects
}
\description{
  Methods for geometrical transformations of
  layered objects (class \code{"layered"}).
}
\usage{
  \method{c}{layered}(\dots)

  \method{shift}{layered}(X, vec=c(0,0), \dots)

  \method{rotate}{layered}(X, \dots, centre=NULL)

  \method{affine}{layered}(X, \dots)

  \method{reflect}{layered}(X)

  \method{flipxy}{layered}(X)

  \method{rescale}{layered}(X, s, unitname)

  \method{scalardilate}{layered}(X, \dots)
}
\arguments{
  \item{X}{
    Object of class \code{"layered"}.
  }
  \item{\dots}{
    Arguments passed to the relevant methods
    when applying the operation to each layer of \code{X}.
  }
  \item{s}{
    Rescaling factor passed to the relevant method for
    \code{\link{rescale}}. May be missing.
  }
  \item{vec}{
    Shift vector (numeric vector of length 2).
  }
  \item{centre}{
    Centre of rotation.
    Either a vector of length 2, or a character string
    (partially matched to \code{"centroid"}, \code{"midpoint"}
    or \code{"bottomleft"}).
    The default is the coordinate origin \code{c(0,0)}.
  }
  \item{unitname}{
    Optional. New name for the unit of length.
    A value acceptable to the function \code{\link{unitname<-}}
  }
}
\details{
  These are methods for the generic functions
  \code{\link{c}},
  \code{\link{shift}},
  \code{\link{rotate}},
  \code{\link{reflect}},
  \code{\link{affine}},
  \code{\link{rescale}},
  \code{\link{scalardilate}} and
  \code{\link{flipxy}}
  for the class of layered objects.

  A layered object represents data that should be plotted in
  successive layers, for example, a background and a foreground.
  See \code{\link{layered}}.

  The method for \code{c} can be used to concatenate two or more
  layered objects \code{\dots} into a single layered object, retaining the plot
  arguments.

  The other methods apply geometrical operations to each of the layers
  of the layered object \code{X}.
}
\value{
  Another object of class \code{"layered"}.
}
\author{
  \adrian
  and \rolf
}
\seealso{
  \code{\link{layered}}
}
\examples{
  B <- owin(c(5500, 9000), c(2500, 7400))
  L <- layered(Window(demopat), unmark(demopat)[B])
  plot(L)
  plot(rotate(L, pi/4))
}
\keyword{spatial}
\keyword{methods}
