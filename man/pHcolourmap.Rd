\name{pHcolourmap}
\alias{pHcolourmap}
\title{
  Colour Map for pH Values
}
\description{
  Create a colour map for values of \eqn{pH}.
}
\usage{
 pHcolourmap(range = c(0, 14), ..., n=256, step = FALSE)
}
\arguments{
  \item{n}{
    Number of different colour values to be used, when \code{step=FALSE}.
  }
  \item{range}{
    Range of \eqn{pH} values that will be accepted as inputs to the
    colour map. A numeric vector of length 2 giving the minimum and
    maximum values of \eqn{pH}.
  }
  \item{step}{
    Logical value. If \code{step=FALSE} (the default) the colours change
    continuously with increasing values of the input. If
    \code{step=TRUE}, the colour is constant on each unit interval of
    \eqn{pH} values.
  }
  \item{\dots}{Ignored.}
}
\details{
  In chemistry the hydrogen potential \eqn{pH} measures how acidic or basic
  a solution is. This function produces colour maps for numerical values
  of \eqn{pH}, always using the same consistent mapping of \eqn{pH}
  values to colours: neutral \eqn{pH=7} is green,
  acidic values \eqn{pH < 7} are yellow or red,
  and basic values \eqn{pH > 7} are blue.

  The argument \code{range} specifies the range of \eqn{pH} values that
  will be mapped by the resulting colour map. It should be a numeric
  vector of length 2 giving the minimum and maximum values of \eqn{pH}
  that the colour map will handle. (Colour maps created with different
  values of \code{range} use essentially the same mapping of colours,
  but when plotted as colour ribbons, display only the specified range.)

  If \code{step=FALSE} (the default) the colours change
  continuously with increasing values of the input. There will be
  \code{n} different colour values in the colour map. Usually \code{n}
  should be a large number.

  If \code{step=TRUE}, the colour is constant on each unit interval of
  \eqn{pH} values. That is, any value of \eqn{pH} in the interval
  \eqn{[k, k+1]}, where \eqn{k} is an integer, will be mapped to the
  same colour.
}
\value{
  A colour map (object of class \code{"colourmap"}).
}

\author{
  \adrian.
}
\seealso{
  \code{\link{colourmap}}
}
\examples{
  plot(pHcolourmap())
  plot(pHcolourmap(step=TRUE))
  plot(pHcolourmap(c(3, 8)))
}
\keyword{spatial}
\keyword{color}
