\name{NAobject}
\alias{NAobject}
\title{
  Make an NA Object
}
\description{
  Create a object which represents a missing or unavailable object
  of a particular class in the \pkg{spatstat} package family.
}
\usage{
  NAobject(cls)
}
\arguments{
  \item{cls}{
    Character string specifying a class.
  }
}
\details{
  This function creates an object which represents
  a missing or unavailable object belonging to one of the
  classes in the \pkg{spatstat} package family.

  For example \code{NAobject("ppp")} creates an object which inherits
  the classes \code{"NAobject"} and \code{"ppp"}, and effectively
  represents a point pattern dataset that is missing or unavailable.

  Typical uses of NA objects include simulation algorithms (where the
  algorithm may sometimes fail to produce a result) and designed
  experiments (where data may be missing for some of the experimental
  cases).

  An NA object can be recognised by \code{\link{is.NAobject}}.

  NA entries in a hyperframe can be recognised by
  \code{\link{is.na.hyperframe}}.
  NA entries in a spatial object list can be recognised by
  \code{\link{is.na.solist}}.
}
\value{
  An object belonging to the classes \code{"NAobject"} and \code{cls}.
}
\author{
  \adrian and \ege.
}
\seealso{
  \code{\link{is.NAobject}},
  \code{\link{is.na.hyperframe}},
  \code{\link{is.na.solist}}.
}
\examples{
  NAobject("tess")
}
\keyword{methods}
