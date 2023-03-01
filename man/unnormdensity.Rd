\name{unnormdensity}
\alias{unnormdensity}
\title{
  Weighted kernel smoother
}
\description{
  An unnormalised version of kernel density estimation
  where the weights are not required to sum to 1.
  The weights may be positive, negative or zero.
}
\usage{
unnormdensity(x, ..., weights = NULL, defaults,
                      weightrange=range(weights))
}
\arguments{
  \item{x}{
    Numeric vector of data
  }
  \item{\dots}{
    Optional arguments passed to \code{\link[stats]{density.default}}.
    Arguments must be \emph{named}.
  }`
  \item{weights}{
    Optional numeric vector of weights for the data.
    The default is equivalent to assuming a weight of 1 for each observation.
  }
  \item{defaults}{
    Optional, named list of arguments passed to
    \code{\link[stats]{density.default}}. These will be overridden by
    arguments in \code{\dots}.
  }
  \item{weightrange}{
    Range of values of \code{weights}, if already known.
  }
}
\details{
  This is an alternative to the standard \R kernel density estimation function
  \code{\link[stats]{density.default}}.

  The standard \code{\link[stats]{density.default}}
  requires the \code{weights} to be nonnegative numbers that add up to 1,
  and returns a probability density (a function that integrates to 1).
  
  This function \code{unnormdensity} does not impose any requirement
  on the \code{weights} except that they be finite. Individual weights may be
  positive, negative or zero. The result is a function that does not
  necessarily integrate to 1 and may be negative. The result is
  the convolution of the kernel \eqn{k} with the weighted data,
  \deqn{
    f(x) = \sum_i w_i k(x- x_i)
  }{
    f(x) = sum of w[i] * k(x - x[i])
  }
  where \eqn{x_i}{x[i]} are the data points and \eqn{w_i}{w[i]} are the
  weights.

  The argument \code{weights} should be a numeric vector of the same
  length as \code{x}, or a single numeric value. The default is to
  assume a weight of 1 for each observation in \code{x}.
  
  The algorithm first selects the kernel bandwidth by
  applying \code{\link[stats]{density.default}} to the data
  \code{x} with normalised, positive weight vector
  \code{w = abs(weights)/sum(abs(weights))} and
  extracting the selected bandwidth.
  Then the result is computed by applying
  applying \code{\link[stats]{density.default}} to \code{x} twice
  using the normalised positive and negative parts of the weights.
  
  Note that the arguments \code{\dots} must be passed by name,
  i.e. in the form (\code{name=value}). Arguments that do not match
  an argument of \code{\link[stats]{density.default}} will be ignored
  \emph{silently}.
}
\section{Warning}{
  If \code{weights} is not specified,
  the default is to assign a weight \eqn{w_i=1}{w[i]=1} to each
  observation \eqn{x_i}{x[i]}.

  This is not the same behaviour as in \code{\link{density.default}} which
  effectively assumes a weight of \eqn{1/n} for each observation \eqn{x_i}{x[i]}
  where \code{n=length(x)}.
}
\value{
  Object of class \code{"density"} as described in
  \code{\link[stats]{density.default}}.
}
\author{
  \adrian
  and \rolf
}
\seealso{
 \code{\link[stats]{density.default}}
}
\examples{
  d <- unnormdensity(1:3, weights=c(-1,0,1), bw=0.3)
  if(interactive()) plot(d)
}
\keyword{smooth}
