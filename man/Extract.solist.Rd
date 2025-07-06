\name{Extract.solist}
\alias{[.solist}
\alias{[<-.solist}
\alias{[[<-.solist}
\title{Extract or Replace Subset of a List of Spatial Objects}
\description{
  Extract or replace some entries in a list of spatial objects,
  or extract a designated sub-region in each object.
}
\usage{
  \method{[}{solist}(x, i, \dots)

  \method{[}{solist}(x, i) <- value

  \method{[[}{solist}(x, i) <- value
}
\arguments{
  \item{x}{
    An object of class \code{"solist"} representing a list of
    two-dimensional spatial objects.
  }
  \item{i}{
    Subset index. Any valid subset index for vectors in the usual \R sense,
    or a window (object of class \code{"owin"}).
  }
  \item{value}{
    Replacement value for the subset. See Details.
  }
  \item{\dots}{Ignored.}
}
\value{
  Another object of the same class as \code{x}.
}
\details{
  These are methods for extracting and replacing subsets
  for the class \code{"solist"}.

  The argument \code{x} should be an object of class \code{"solist"}
  representing a list of two-dimensional spatial objects.
  See \code{\link{solist}}.

  For the subset method, the subset index \code{i} can be either
  a vector index (specifying some elements of the list)
  or a spatial window (specifying a spatial sub-region). 

  For the replacement \code{"[<-"} method,
  \code{i} must be a vector index and \code{value} must be a list
  of spatial objects;
  the designated elements of \code{x} will be
  replaced by the corresponding elements of the list \code{value}.

  For the replacement \code{"[[<-"} method,
  \code{i} must be a single integer, and \code{value} must be a spatial
  object.

  If all entries of \code{x} belong to the same class,
  then in the replacement methods, a value of \code{NA}
  is automatically coerced to an NA object of the same class.
}
\seealso{
  \code{\link{solist}},
  \code{\link{plot.solist}},
  \code{\link{summary.solist}}
}
\examples{
   x <- solist(japanesepines, cells, redwood)
   x[2:3]
   x[square(0.5)]
   x[1] <- list(finpines)
   x[[1]] <- cells
   x[[3]] <- NA
 }
\author{
  \spatstatAuthors
}
\keyword{spatial}
\keyword{list}
\keyword{manip}
