\name{pairdist}
\alias{pairdist}
\title{Pairwise distances}
\description{
  Computes the matrix of distances between all pairs of `things'
  in a dataset
}
\usage{
  pairdist(X, \dots)
}
\arguments{
  \item{X}{
    Object specifying the locations of a set of `things'
    (such as a set of points or a set of line segments).
  }
  \item{\dots}{
    Further arguments depending on the method.
  }
}
\value{
  A square matrix whose \code{[i,j]} entry is the distance
  between the `things' numbered \code{i} and \code{j}.
}
\details{
  Given a dataset \code{X} and \code{Y}
  (representing either a point pattern or
  a line segment pattern)
  \code{pairdist} computes the distance between each pair of
  `things' in the dataset, and returns 
  a matrix containing these distances.

  The function \code{pairdist} is generic, with
  methods for point patterns (objects of class \code{"ppp"}),
  line segment patterns (objects of class \code{"psp"})
  and a default method. See the documentation
  for \code{\link{pairdist.ppp}}, \code{\link{pairdist.psp}}
  or \code{\link{pairdist.default}} for details.
}
\section{Distance values}{
  The values returned by \code{pairdist(X)} are distances,
  expressed as multiples of the unit of length of the spatial
  coordinates in \code{X}.
  The unit of length is given by \code{\link{unitname}(X)}.

  Note that, if the unit of length in \code{X} is a composite
  expression such as \sQuote{2 microns},
  then the values of \code{pairdist(X)} are expressed as multiples
  of 2 microns, rather than being expressed in microns.
}
\seealso{
  \code{\link{pairdist.ppp}},
  \code{\link{pairdist.psp}},
  \code{\link{pairdist.default}},
  \code{\link{crossdist}},
  \code{\link{nndist}},
  \code{\link[spatstat.explore]{Kest}}
}
\author{
  \pavel
  and 
  \adrian.
}
\keyword{spatial}
\keyword{math}

