\name{weighted.median}
\alias{weighted.median}
\alias{weighted.quantile}
\alias{weighted.var}
\title{
  Weighted Median, Quantiles or Variance
}
\description{
  Compute the median, quantiles or variance of a set of numbers which have
  weights associated with them.
}
\usage{
weighted.median(x, w, na.rm = TRUE, type=2)

weighted.quantile(x, w, probs=seq(0,1,0.25), na.rm = TRUE, type=4)

weighted.var(x, w, na.rm = TRUE)
}
\arguments{
  \item{x}{
    Data values.
    A vector of numeric values, for which the median or quantiles are required.
  }
  \item{w}{
    Weights.
    A vector of nonnegative numbers, of the same length as \code{x}.
  }
  \item{probs}{
    Probabilities for which the quantiles should be computed.
    A numeric vector of values between 0 and 1.
  }
  \item{na.rm}{
    Logical. Whether to ignore \code{NA} values.
  }
  \item{type}{
    Integer specifying the rule for calculating the median or quantile,
    corresponding to the rules available for
    \code{\link[stats]{quantile}}.
    The only valid choices are \code{type=1}, \code{2} or \code{4}.
    See Details.
  }
    
}
\details{
  The \code{i}th observation \code{x[i]} is treated as having
  a weight proportional to \code{w[i]}.

  The weighted median is a value \code{m}
  such that the total weight of data less than or equal to \code{m}
  is equal to half the total weight. More generally, the weighted quantile with
  probability \code{p} is a value \code{q} 
  such that the total weight of data less than or equal to \code{q}
  is equal to \code{p} times the total weight.
  
  If there is no such value, then
  \itemize{
    \item if \code{type=1}, the next largest value is returned
    (this is the right-continuous inverse of the left-continuous
    cumulative distribution function);
    \item if \code{type=2}, the average of the two surrounding values is
    returned (the average of the right-continuous and left-continuous
    inverses);
    \item if \code{type=4}, linear interpolation is performed.
  }

  Note that the default rule for \code{weighted.median} is
  \code{type=2}, consistent with the traditional definition of the median, 
  while the default for \code{weighted.quantile} is \code{type=4}.
}
\value{
  A numeric value or vector.
}
\author{
  \adrian.
}
\seealso{
  \code{\link[stats]{quantile}}, \code{\link[stats]{median}}.
}
\examples{
  x <- 1:20
  w <- runif(20)
  weighted.median(x, w)
  weighted.quantile(x, w)
  weighted.var(x, w)
}
\keyword{math}
