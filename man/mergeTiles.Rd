\name{mergeTiles}
\alias{mergeTiles}
\alias{mergeTiles.tess}
\title{
  Merge Groups of Tiles in a Tessellation
}
\description{
  Given a tessellation and a factor which classifies the tiles of the
  tessellation into groups, form the set union of each group of tiles,
  and construct the tessellation based on these merged tiles.
}
\usage{
 mergeTiles(x, group, \dots)

 \method{mergeTiles}{tess}(x, group, \dots)
}
\arguments{
  \item{x}{
    A tessellation (object of class \code{"tess"}).
  }
  \item{group}{
    A \code{factor}, with one entry for each tile of \code{x},
    which classifies the tiles into groups.
  }
  \item{\dots}{Other arguments, currently ignored.}
}
\value{
  A tessellation.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{tess}}
}
\examples{
  D <- dirichlet(cells)
  g <- factor(sample(letters[1:4], npoints(cells), replace=TRUE))
  plot(D, values=g, main="tiles to be merged")
  plot(mergeTiles(D, g), do.col=TRUE, "merged")
}
\keyword{spatial}
\keyword{manip}
