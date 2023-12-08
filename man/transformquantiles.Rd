\name{transformquantiles}
\alias{transformquantiles}
\title{
 Transform the Quantiles
}
\description{
  Apply a transformation to the quantiles of a vector,
  or to the quantiles of the pixel values in a pixel image.
}
\usage{
transformquantiles(X, uniform = FALSE, reverse = FALSE, ...)
}
\arguments{
  \item{X}{
    A numeric vector, matrix, array, or a pixel image
    (object of class \code{"im"}).
  }
  \item{uniform}{
    Logical value specifying whether each quantile value should be
    replaced by the corresponding cumulative probability
    (called \emph{histogram equalisation}, \emph{transformation to uniformity}
    or \emph{probability integral transformation}).
  }
  \item{reverse}{
    Logical value specifying whether to swap the
    upper and lower quantiles. 
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  The argument \code{X} may be a vector, matrix,
  array, or a pixel image (object of class \code{"im"}).

  The algorithm will first extract the entries or pixel values of
  \code{X} as a vector, and sort the values into ascending order.

  If \code{uniform=TRUE}, the entries in this vector will be replaced by the
  corresponding cumulative probabilities (the \code{k}th
  smallest value will be replaced by the
  number \code{(k-0.5)/n} where \code{n} is the total number of values).

  If \code{reverse=TRUE}, the resulting vector will be reversed
  so that it is in descending order (so that the \code{k}th smallest
  value will be swapped with the \code{k}th largest value).

  Finally the transformed values will be replaced into the original
  positions in the vector, matrix, array, or pixel image.

  The case \code{uniform=TRUE, reverse=FALSE} is called
  \emph{transformation to uniformity}, the
  \emph{probability integral transformation},
  \emph{histogram equalisation}, or \emph{quantile transformation}.
  The resulting values are uniformly distributed between 0 and 1;
  a histogram of the values in \code{X} is flat.
}
\value{
  Another object of the same type as \code{X}.
}
\author{
  \spatstatAuthors.
}
\seealso{
  To apply an arbitrary function \code{f} to the pixel values in an image,
  use the idiom \code{X[] <- f(X[])}.
%  or use \code{\link{eval.im}}.
}
\examples{
  X <- c(3, 5, 1, 2, 4)
  transformquantiles(X, reverse=TRUE)
  transformquantiles(X, uniform=TRUE)
  transformquantiles(X, uniform=TRUE, reverse=TRUE)

%  Z <- 100 * as.im(bei.extra$grad, dimyx=c(25,50))
%  opa <- par(mfrow=c(2,2))
%  plot(Z,main="original")
%  plot(transformquantiles(Z, reverse=TRUE), main="reverse=TRUE")
%  plot(transformquantiles(Z, uniform=TRUE),
%  main="uniform=TRUE", col=grey.colors,zlim=c(0,1))
%  plot(transformquantiles(Z, uniform=TRUE, reverse=TRUE),
%                          main="uniform=TRUE, reverse=TRUE",
%                          col=grey.colors, zlim=c(0,1))
%  par(opa)
}
\keyword{spatial}
\keyword{manip}
\keyword{nonparametric}
\keyword{univar}
\concept{Histogram equalisation}
\concept{Quantile transformation}
\concept{Probability integral transformation}
