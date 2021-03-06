\name{edges}
\alias{edges}
\title{
  Extract Boundary Edges of a Window.
}
\description{
  Extracts the boundary edges of a window
  and returns them as a line segment pattern.
}
\usage{
  edges(x, \dots, window = NULL, check = FALSE) 
}
\arguments{
  \item{x}{
    A window (object of class \code{"owin"}), or
    data acceptable to \code{\link{as.owin}}, specifying the window
    whose boundary is to be extracted.
  }
  \item{\dots}{
    Ignored.
  }
  \item{window}{
    Window to contain the resulting line segments.
    Defaults to \code{as.rectangle(x)}.
  }
  \item{check}{
    Logical. Whether to check the validity of the resulting
    segment pattern.
  }
}
\details{
  The boundary edges of the window \code{x} will be extracted as a line
  segment pattern.
}
\value{
  A line segment pattern (object of class \code{"psp"}).
}
\seealso{
  \code{\link{perimeter}} for calculating the total length of the boundary.
}
\examples{
  edges(square(1))
  edges(letterR)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{manip}
