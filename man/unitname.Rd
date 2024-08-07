\name{unitname}  
\alias{unitname}
\alias{unitname.im}
\alias{unitname.owin}
\alias{unitname.ppp}
\alias{unitname.psp}
\alias{unitname.quad}
\alias{unitname.tess}
\alias{unitname<-}
\alias{unitname<-.im}
\alias{unitname<-.owin}
\alias{unitname<-.ppp}
\alias{unitname<-.psp}
\alias{unitname<-.quad}
\alias{unitname<-.tess}
\title{Name for Unit of Length}
\description{
  Inspect or change the name of the unit of length
  in a spatial dataset.
}
\usage{
unitname(x)
unitname(x) <- value
\method{unitname}{im}(x)
\method{unitname}{owin}(x)
\method{unitname}{ppp}(x)
\method{unitname}{psp}(x)
\method{unitname}{quad}(x)
\method{unitname}{tess}(x)
\method{unitname}{im}(x) <- value
\method{unitname}{owin}(x) <- value
\method{unitname}{ppp}(x) <- value
\method{unitname}{psp}(x) <- value
\method{unitname}{quad}(x) <- value
\method{unitname}{tess}(x) <- value
}
\arguments{
  \item{x}{A spatial dataset.
    Either a point pattern (object of class \code{"ppp"}),
    a line segment pattern (object of class \code{"psp"}),
    a window (object of class \code{"owin"}),
    a pixel image (object of class \code{"im"}),
    a tessellation (object of class \code{"tess"}),
    a quadrature scheme (object of class \code{"quad"}),
    or a fitted point process model
    (object of class \code{"ppm"} or \code{"kppm"} or \code{"slrm"}
    or \code{"dppm"} or \code{"minconfit"}).
  }
  \item{value}{
    Name of the unit of length. See Details.
  }
}
\details{
  Spatial datasets in the \pkg{spatstat} package
  may include the name of the unit of length. This name is used
  when printing or plotting the dataset, and in some other
  applications. 

  \code{unitname(x)} extracts this name,
  and \code{unitname(x) <- value} sets the name to \code{value}.

  A valid name is either
  \itemize{
    \item a single character string
    \item a vector of two character strings giving the
    singular and plural forms of the unit name
    \item a list of length 3, containing two character strings
    giving the singular and plural forms of the basic unit,
    and a number specifying the multiple of this unit.
  }
  
 Note that re-setting the name of the unit of length \emph{does not}
 affect the numerical values in \code{x}. It changes only the string
 containing the name of the unit of length. To rescale the numerical
 values, use \code{\link[spatstat.geom]{rescale}}.
}
\value{
  The return value of \code{unitname} is an object of class \code{"unitname"}
  containing the name of the unit of length in \code{x}. There are
  methods for \code{print}, \code{summary}, \code{as.character},
  \code{\link[spatstat.geom]{rescale}} and \code{\link[spatstat.geom]{compatible}}.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link[spatstat.geom]{rescale}},
  \code{\link[spatstat.geom]{owin}},
  \code{\link[spatstat.geom]{ppp}}
}

\examples{
  X <- runifrect(20)

  # if the unit of length is 1 metre:
  unitname(X) <- c("metre", "metres")

  # if the unit of length is 6 inches:
  unitname(X) <- list("inch", "inches", 6)
}
\keyword{spatial}
\keyword{manip}

