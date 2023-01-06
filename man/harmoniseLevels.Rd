\name{harmoniseLevels}
\alias{harmoniseLevels}
\title{
  Harmonise the levels of several factors, or factor-valued pixel images.
}
\description{
  Given several factors (or factor-valued pixel images)
  convert them so that they all use the same set of levels.
}
\usage{
harmoniseLevels(...)
}
\arguments{
  \item{\dots}{
    Factors, or factor-valued pixel images.
  }
}
\details{
  All of the arguments \code{\dots} must be factors, or factor-valued
  pixel images (objects of class \code{"im"}).

  The \code{\link{levels}} of each factor will be extracted, and
  combined by taking the union of all the levels. Then each factor will
  be converted to a new factor so that all of the new factors have
  exactly the same set of levels.
}
\value{
  A list, containing the same number of arguments as the input,
  consisting of factors or factor-valued pixel images.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{levels}}, \code{\link{levels.im}}
}
\examples{
  (a <- factor(sample(letters[1:3], 10, replace=TRUE)))
  (b <- factor(sample(LETTERS[1:4], 7, replace=TRUE)))
  harmoniseLevels(a,b)

  (A <- gorillas.extra$vegetation)
  (B <- gorillas.extra$slopetype)
  harmoniseLevels(A,B)
}
\keyword{manip}
