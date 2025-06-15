\name{nnmark}
\alias{nnmark}
\title{
  Mark of Nearest Neighbour 
}
\description{
  Given a marked point pattern dataset \code{X}
  this function computes, for each desired location \code{y},
  the mark attached to the nearest neighbour of \code{y} in \code{X}.
  The desired locations \code{y} can be either a pixel grid
  or the point pattern \code{X} itself.
}
\usage{
nnmark(X, \dots, k = 1, at=c("pixels", "points"),
                 ties=c("first", "mean", "min", "max"),
                 proper=FALSE)
}
\arguments{
  \item{X}{
    A marked point pattern (object of class \code{"ppp"}).
  }
  \item{\dots}{
    Arguments passed to \code{\link[spatstat.geom]{as.mask}} to determine the
    pixel resolution.
  }
  \item{k}{
    Single integer. The \code{k}th nearest data point will be used.
  }
  \item{at}{
    String specifying whether to compute the values
    at a grid of pixel locations (\code{at="pixels"}) or
    only at the points of \code{X} (\code{at="points"}).
  }
  \item{ties}{
    Character string (partially matched) indicating how to handle
    the case of ties, where there are two or more data points at the same
    location. See Details.
  }
  \item{proper}{
    Logical value specifying how to define nearest neighbours
    if there are two or more data points at the same location.
    Applies only when \code{at="points"}.
    See Details.
  }
}
\details{
  Given a marked point pattern dataset \code{X}
  this function computes, for each desired location \code{y},
  the mark attached to the point of \code{X} that is nearest
  to \code{y}. The desired locations \code{y} can be either a pixel grid
  or the point pattern \code{X} itself.

  The argument \code{X} must be a marked point pattern (object
  of class \code{"ppp"}, see \code{\link{ppp.object}}).
  The marks are allowed to be a vector or a data frame.

  \itemize{
    \item
    If \code{at="points"}, then for each point in \code{X},
    the algorithm finds the nearest \emph{other} point in \code{X},
    and extracts the mark attached to it.
    The result is a vector or data frame containing the marks
    of the neighbours of each point.
    \item
    If \code{at="pixels"} (the default), then for each pixel
    in a rectangular grid, the algorithm finds the nearest point in \code{X},
    and extracts the mark attached to it.
    The result is an image or a list of images containing the marks
    of the neighbours of each pixel.
    The pixel resolution is controlled by the arguments \code{\dots}
    passed to \code{\link[spatstat.geom]{as.mask}}.
  }

  If the argument \code{k} is given, then the \code{k}-th nearest
  neighbour will be used.

  The arguments \code{ties} and \code{proper} specify how to handle
  the case where two or more data points are at the same spatial
  location.
  \itemize{
    \item
    \code{ties} determines how to pool the mark values.
    If \code{ties="first"} (the default), the mark value for
    this location is taken to be the mark of the data point that is listed
    first in sequence in the dataset \code{X}. If \code{ties="mean"},
    \code{ties="max"} or \code{ties="min"}, the mark value for this
    location is taken to be the mean, maximum or minimum (respectively) of
    the mark values of all the data points at this location (after
    converting the mark values to numerical values). 
    \item
    \code{proper} determines how to define nearest neighbours,
    when \code{at="points"}.
    If \code{proper=TRUE}, the nearest neighbour of a data point
    must be another data point lying a nonzero distance away from it.
    If \code{proper=FALSE} (the default), then two data points
    occupying the exact same spatial location can be nearest neighbours.
  }
}  
\value{
  \emph{If \code{X} has a single column of marks:}
  \itemize{
    \item 
    If \code{at="pixels"} (the default), the result is
    a pixel image (object of class \code{"im"}). 
    The value at each pixel is the mark attached
    to the nearest point of \code{X}.
    \item
    If \code{at="points"}, the result is a vector or factor
    of length equal to the number of points in \code{X}.
    Entries are the mark values of the
    nearest neighbours of each point of \code{X}.
  }
  \emph{If \code{X} has a data frame of marks:}
  \itemize{
    \item 
    If \code{at="pixels"} (the default), the result is a named list of 
    pixel images (object of class \code{"im"}). There is one
    image for each column of marks. This list also belongs to
    the class \code{"solist"}, for which there is a plot method.
    \item
    If \code{at="points"}, the result is a data frame
    with one row for each point of \code{X},
    Entries are the mark values of the
    nearest neighbours of each point of \code{X}.
  }
}
\author{\spatstatAuthors.}
\seealso{
  \code{\link[spatstat.explore]{Smooth.ppp}},
  \code{\link[spatstat.explore]{marktable}},
  \code{\link{nnwhich}}
}
\examples{
  plot(nnmark(ants))
  v <- nnmark(ants, at="points")
  v[1:10]
  plot(nnmark(finpines))
  vf <- nnmark(finpines, at="points")
  vf[1:5,]
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
