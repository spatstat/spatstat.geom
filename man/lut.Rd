\name{lut}
\alias{lut}
\title{Lookup Tables}
\description{
  Create a lookup table.
}
\usage{
lut(outputs, ..., range=NULL, breaks=NULL, inputs=NULL,
    gamma=1, compress=NULL, decompress=NULL)
}
\arguments{
  \item{outputs}{Vector of output values.}
  \item{\dots}{Ignored.}
  \item{range}{
    Interval of numbers to be mapped.
    A numeric vector of length 2, specifying the ends of the range of values
    to be mapped.
    Incompatible with \code{breaks} or \code{inputs}.
  }
  \item{inputs}{
    Input values to which the output values are associated.
    A factor or vector of the same length as \code{outputs}.
    Incompatible with \code{breaks} or \code{range}.
  }
  \item{breaks}{
    Breakpoints for the lookup table.
    A numeric vector of length equal to \code{length(outputs)+1}.
    Incompatible with \code{range} or \code{inputs}.
  }
  \item{gamma}{
    Exponent for gamma correction, when \code{range} is given.
    A single positive number.
    See Details.
  }
  \item{compress}{
    Optional. Experimental.
    An \R function determining a nonlinear transformation of the
    domain of the lookup table.
    See section on \emph{Nonlinear lookup tables}.
  }
  \item{decompress}{
    Experimental.
    An \R function giving the inverse function of \code{compress},
    if \code{compress} is specified.
  }
}
\details{
  A lookup table is a function, mapping input values to output values.

  The command \code{lut} creates an object representing
  a lookup table, which can then be used to control various behaviour
  in the \pkg{spatstat} package. It can also be used to compute the
  output value assigned to any input value. 

  The argument \code{outputs} specifies the output values to which
  input data values will be mapped. It should be a vector of
  any atomic type (e.g. numeric, logical, character, complex) or factor
  values. 

  Exactly one of the arguments \code{range}, \code{inputs} or \code{breaks}
  must be specified by name.
  \itemize{
    \item 
    If \code{inputs} is given, then it should be a vector or factor,
    of the same length as \code{outputs}. The entries of \code{inputs} can be
    any atomic type (e.g. numeric, logical, character, complex) or factor
    values. The resulting lookup table associates the value \code{inputs[i]}
    with the value \code{outputs[i]}.
    The argument \code{outputs} should have the same length as \code{inputs}.
    \item
    If \code{range} is given, then it determines the interval of the real
    number line that will be mapped. It should be a numeric vector of
    length 2. 
    The interval will be divided evenly into bands, each of which is
    mapped to an entry of \code{outputs}. (If \code{gamma} is given,
    then the bands are equally spaced on a scale where the original values
    are raised to the power \code{gamma}.)
    (See the section on \emph{Nonlinear lookup tables} for the case
    where \code{compress} and \code{decompress} are given.)
    \item
    If \code{breaks} is given, then it determines intervals
    of the real number line
    which are mapped to each output value. It should be a numeric vector,
    of length at least 2, with entries that are in increasing order.
    Infinite values are allowed. Any number in the range
    between \code{breaks[i]} and \code{breaks[i+1]} will be mapped to the
    value \code{outputs[i]}. 
    The argument \code{outputs} should have length equal to
    \code{length(breaks) - 1}.
  }
  It is also permissible for \code{outputs} to be a single value,
  representing a trivial lookup table in which all data values are mapped to
  the same output value.

  The result is an object of class \code{"lut"}. 
  There is a \code{print} method for this class.
  Some plot commands in the \pkg{spatstat} package accept an object
  of this class as a specification of a lookup table.

  The result is also a function \code{f} which can be used to compute
  the output value assigned to any input data value. 
  That is, \code{f(x)} returns the output value assigned
  to \code{x}. This also works for vectors of input data values.
}
\section{Nonlinear lookup tables}{
  If the arguments \code{compress} and \code{decompress} are given,
  they define a transformation of the range of numbers.
  A typical example would be a logarithmic lookup table defined
  by \code{compress = log10} and
  \code{decompress = function(x) { 10^xg }}.

  These functions have no effect on the interpretation of the
  arguments \code{range}, \code{breaks} and \code{inputs}. However,
  if \code{range} is given, then the range of values will be divided
  into intervals which have equal length on the scale defined by
  \code{compress}. That is, the range of numbers determined by
  \code{compress(range)} will be evenly divided into intervals,
  and these intervals will be mapped back to the original scale by
  \code{decompress} to determine the \code{breaks}.
  For a logarithmic lookup table with \code{compress=log10},
  the specified range will be divided into intervals which are equal
  on a logarithmic scale.

  The arguments \code{compress} and \code{decompress}
  should be functions which are vectorised (i.e. if \code{x} is a
  vector then \code{compress(x)} and
  \code{decompress(x)} are also vectors of the same length as \code{x})
  and increasing (if \code{x < y} then \code{compress(x) < compress(y)}
  and \code{decompress(x) < decompress(y)}.
}
\value{
  A function, which is also an object of class \code{"lut"}.
}
\seealso{
  \code{\link{colourmap}}.
}
\examples{
  ## lookup table for real numbers, using breakpoints
  cr <- lut(factor(c("low", "medium", "high")), breaks=c(0,5,10,15))
  cr
  cr(3.2)
  cr(c(3,5,7))
  ## lookup table for discrete set of values
  ct <- lut(c(0,1), inputs=c(FALSE, TRUE))
  ct(TRUE)
  ## logarithmic lookup table
  cl <- lut(letters[1:3], range=c(0.1, 100),
            compress=log10, decompress=function(x) { 10^x })
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{manip}
