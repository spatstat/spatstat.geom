\name{layout.boxes}
\alias{layout.boxes}
\title{
  Generate a Row or Column Arrangement of Rectangles.
}
\description{
  A simple utility to generate a row or column of boxes (rectangles)
  for use in point-and-click panels.
}
\usage{
layout.boxes(B, n, horizontal = FALSE, aspect = 0.5, usefrac = 0.9)
}
\arguments{
  \item{B}{
    Bounding rectangle for the boxes.
    An object of class \code{"owin"}.
  }
  \item{n}{
    Integer. The number of boxes.
  }
  \item{horizontal}{
    Logical. If \code{TRUE}, arrange the boxes in a horizontal row.
    If \code{FALSE} (the default), arrange them in a vertical column.
  }
  \item{aspect}{
    A single finite positive number, giving the
    aspect ratio (height divided by width) of each box,
    or \code{NA} or \code{Inf},
    indicating that the aspect ratio is unconstrained.
  }
  \item{usefrac}{
    Number between 0 and 1. The 
    fraction of height or width of \code{B} that should be
    occupied by boxes.
  }
}
\details{
  This simple utility generates a list of boxes (rectangles)
  inside the bounding box \code{B} arranged in a regular
  row or column. It is useful for generating the
  positions of the panel buttons in the function
  \code{\link{simplepanel}}.

  The argument \code{aspect} specifies the ratio of height to width
  (height divided by width). 
  If \code{aspect} is a finite numerical value, then the boxes will
  have the given aspect ratio.
  If \code{aspect} is \code{Inf} or
  \code{NA}, aspect ratio is unconstrained; the boxes will have the
  maximum possible width and height.
}
\value{
  A list of rectangles (objects of class \code{"owin"} which are rectangles).
}
\examples{
  B <- owin(c(0,10),c(0,1))
  boxes <- layout.boxes(B, 5, horizontal=TRUE)
  plot(B, main="", col="blue")
  niets <- lapply(boxes, plot, add=TRUE, col="grey")
}
\author{
  \adrian
  and \rolf
}
\seealso{
  \code{\link{simplepanel}}
}
\keyword{utilities}


