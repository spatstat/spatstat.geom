\name{colourmap}
\alias{colourmap}
\title{Colour Lookup Tables}
\description{
  Create a colour map (colour lookup table).
}
\usage{
colourmap(col, \dots, range=NULL, breaks=NULL, inputs=NULL, gamma=1,
          compress=NULL, decompress=NULL)
}
\arguments{
  \item{col}{Vector of values specifying colours}
  \item{\dots}{Ignored.}
  \item{range}{
    Interval to be mapped.
    A numeric vector of length 2, specifying the endpoints of the
    range of values to be mapped.
    Incompatible with \code{breaks} or \code{inputs}.
  }
  \item{inputs}{
    Values to which the colours are associated.
    A factor or vector of the same length as \code{col}.
    Incompatible with \code{breaks} or \code{range}.
  }
  \item{breaks}{
    Breakpoints for the colour map.
    A numeric vector of length equal to \code{length(col)+1}.
    Incompatible with \code{range} or \code{inputs}.
  }
  \item{gamma}{
    Exponent for the gamma correction, when \code{range} is given.
    A single positive number.
    See Details.
  }
  \item{compress}{
    Optional. Experimental.
    An \R function determining a nonlinear transformation of the
    domain of the colour map.
    See section on \emph{Nonlinear colour maps}.
  }
  \item{decompress}{
    Experimental.
    An \R function giving the inverse function of \code{compress},
    if \code{compress} is specified.
  }
}
\details{
  A colour map is a mechanism for associating colours with data.
  It can be regarded as a function, mapping data to colours.

  The command \code{colourmap} creates an object representing
  a colour map, which can then be used to control the plot commands
  in the \pkg{spatstat} package. It can also be used to compute the
  colour assigned to any data value. 

  The argument \code{col} specifies the colours to which
  data values will be mapped. It should be a vector
  whose entries can be interpreted as colours by the standard
  \R graphics system. The entries can be string names of colours
  like \code{"red"}, or integers that refer to
  colours in the standard palette, or strings containing
  six-letter hexadecimal codes like \code{"#F0A0FF"}. 

  Exactly one of the arguments \code{range}, \code{inputs} or \code{breaks}
  must be specified by name.

  \itemize{
    \item
    If \code{inputs} is given, then it should be a vector or factor,
    of the same length as \code{col}. The entries of \code{inputs} can be
    any atomic type (e.g. numeric, logical, character, complex) or factor
    values. The resulting colour map associates the value \code{inputs[i]}
    with the colour \code{col[i]}.
    The argument \code{col} should have the same length as \code{inputs}.
    \item
    If \code{range} is given, then it determines the interval of the real
    number line that will be mapped. It should be a numeric vector of
    length 2. The interval will be divided evenly into bands, each of which is
    assigned one of the colours in \code{col}. (If \code{gamma} is given,
    then the bands are equally spaced on a scale where the original values
    are raised to the power \code{gamma}.)
    (See the section on \emph{Nonlinear colour maps} for the case
    where \code{compress} and \code{decompress} are given.)
    \item
    If \code{breaks} is given, then it determines the precise intervals
    of the real number line
    which are mapped to each colour. It should be a numeric vector,
    of length at least 2, with entries that are in increasing order.
    Infinite values are allowed. Any number in the range
    between \code{breaks[i]} and \code{breaks[i+1]} will be mapped to the
    colour \code{col[i]}. 
    The argument \code{col} should have length equal to
    \code{length(breaks) - 1}.
  }
  It is also permissible for \code{col} to be a single colour value,
  representing a trivial colour map in which all data values are mapped to
  the same colour.
  
  The result is an object of class \code{"colourmap"}. 
  There are \code{print} and \code{plot} methods for this class.
  Some plot commands in the \pkg{spatstat} package accept an object
  of this class as a specification of the colour map.

  The result is also a function \code{f} which can be used to compute
  the colour assigned to any data value. 
  That is, \code{f(x)} returns the character value of the colour assigned
  to \code{x}. This also works for vectors of data values.
}
\section{Nonlinear colour maps}{
  If the arguments \code{compress} and \code{decompress} are given,
  they define a transformation of the range of numbers.
  A typical example would be a logarithmic colour map defined
  by \code{compress = log10} and
  \code{decompress = function(x) { 10^x }}.

  These functions have no effect on the interpretation of the
  arguments \code{range}, \code{breaks} and \code{inputs}. However,
  if \code{range} is given, then the range of values will be divided
  into intervals which have equal length on the scale defined by
  \code{compress}. That is, the range of numbers determined by
  \code{compress(range)} will be evenly divided into intervals,
  and these intervals will be mapped back to the original scale by
  \code{decompress} to determine the \code{breaks}.
  For a logarithmic colour map with \code{compress=log10},
  the specified range will be divided into intervals which are equal
  on a logarithmic scale.

  The arguments \code{compress} and \code{decompress} affect the way
  in which the colour map is plotted by \code{\link{plot.colourmap}}.
  For a continuous colour map, the range of input values is plotted on
  the compressed scale, but annotated on the original scale.
  See the Examples.

  The arguments \code{compress} and \code{decompress}
  should be functions which are vectorised (i.e. if \code{x} is a
  vector then \code{compress(x)} and
  \code{decompress(x)} are also vectors of the same length as \code{x})
  and increasing (if \code{x < y} then \code{compress(x) < compress(y)}
  and \code{decompress(x) < decompress(y)}.
}
\value{
  A function, which is also an object of class \code{"colourmap"}
  and \code{"lut"}.
}
\seealso{
  The plot method \code{\link{plot.colourmap}}.
  
  See the \R help file on
  \code{\link[grDevices:colors]{colours}} for information about the colours
  that \R recognises, and how to manipulate them.
  
  To make a smooth transition between colours, see
  \code{\link{interp.colourmap}}.
  To alter individual colour values, see
  \code{\link{tweak.colourmap}}.
  To extract or replace all colour values, see
  \code{\link{colouroutputs}}.

  See also \code{\link{restrict.colourmap}}
  and \code{\link{rev.colourmap}}.
  
  See \code{\link[spatstat.geom:colourtools]{colourtools}}
  for more tools to manipulate colour values.

  See \code{\link{lut}} for lookup tables.
}
\examples{
  ## colour map for real numbers, using breakpoints
  cr <- colourmap(c("red", "blue", "green"), breaks=c(0,5,10,15))
  cr
  cr(3.2)
  cr(c(3,5,7))
  ## a large colour map
  co <- colourmap(rainbow(100), range=c(-1,1))
  co(0.2)
  ## colour map for discrete set of values
  ct <- colourmap(c("red", "green"), inputs=c(FALSE, TRUE))
  ct(TRUE)
  ## logarithmic colour map
  cl <- colourmap(rainbow(25), range=c(0.1, 1000),
                  compress=log10, decompress=function(x) { 10^x })
  plot(cl)
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{color}
