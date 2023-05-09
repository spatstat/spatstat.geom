\name{methods.boxx}
\Rdversion{1.1}
\alias{methods.boxx}  %DoNotExport
\alias{print.boxx}
\alias{unitname.boxx}
\alias{unitname<-.boxx}
\alias{scale.boxx}
\title{
  Methods for Multi-Dimensional Box
}
\description{
  Methods for class \code{"boxx"}.
}
\usage{
  \method{print}{boxx}(x, ...)
  \method{unitname}{boxx}(x)
  \method{unitname}{boxx}(x) <- value
  \method{scale}{boxx}(x, center=TRUE, scale=TRUE)
}
\arguments{
  \item{x}{
    Object of class \code{"boxx"} representing a multi-dimensional box.
  }
  \item{\dots}{
    Other arguments passed to \code{print.default}.
  }
  \item{value}{
    Name of the unit of length. See \code{\link{unitname}}.
  }
  \item{center,scale}{
    Arguments passed to \code{\link[base]{scale.default}}
    to determine the rescaling.
  }
}
\details{
  These are methods for the generic functions
  \code{\link{print}}, \code{\link{unitname}},
  \code{\link{unitname<-}} and \code{\link[base]{scale}}
  for the class \code{"boxx"} of multi-dimensional boxes.

  The 
  \code{print} method prints a description of the box,
  the \code{unitname} method extracts the name of the unit of
  length in which the box coordinates are expressed, while the assignment method
  for \code{unitname} assigns this unit name.

  The \code{scale} method rescales each spatial coordinate of \code{x}.
}
\value{
  For \code{print.boxx} the value is \code{NULL}.
  For \code{unitname.boxx} an object of class \code{"units"}.
  For \code{unitname<-.boxx} and
  \code{scale.boxx} the result is the updated \code{"boxx"} object \code{x}.
}
\author{
  \spatstatAuthors
}
\seealso{
  \code{\link{boxx}},
  \code{\link{is.boxx}},
  \code{\link{print}},
  \code{\link{unitname}},
  \code{\link[base]{scale}}
}
\examples{
   B <- boxx(c(0,10),c(0,10),c(0,5),c(0,1), unitname=c("metre", "metres"))
   B
   is.boxx(B)
   unitname(B)
   # Northern European usage
   unitname(B) <- "meter"
   scale(B)
}
\keyword{spatial}
\keyword{methods}
