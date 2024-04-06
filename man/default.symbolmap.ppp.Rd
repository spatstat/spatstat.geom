\name{default.symbolmap.ppp}
\alias{default.symbolmap.ppp}
\title{
  Default Symbol Map for Point Pattern
}
\description{
  Determines a symbol map for plotting the spatial point pattern \code{x}.
}
\usage{
\method{default.symbolmap}{ppp}(x, \dots,
      chars = NULL, cols = NULL, 
      fixsize = FALSE,
      maxsize = NULL, meansize = NULL, markscale = NULL,
      minsize = NULL, zerosize = NULL, marktransform = NULL)
}
\arguments{
  \item{x}{
    A spatial point pattern (object of class \code{"ppp"}).
  }
  \item{\dots}{
    extra graphical parameters, passed to \code{\link{symbolmap}}
    (and ultimately to \code{\link{points}} and/or
    \code{\link{symbols}}).
  }
  \item{chars}{
    the plotting character(s) used to plot points.
    Either a single character, an integer,
    or a vector of single characters or integers.
    Ignored if \code{symap} is given.
  }
  \item{cols}{
    the colour(s) used to plot points.
    Either an integer index from 1 to 8 (indexing the standard colour
    palette), a character string giving the name of a colour,
    or a string giving the hexadecimal representation of a colour,
    or a vector of such integers or strings.
    See the section on \emph{Colour Specification}
    in the help for \code{\link[graphics]{par}}.
  }
  \item{fixsize}{
    Logical value specifying whether the symbols should all have the
    same physical size on the plot. Default is \code{FALSE}.
  }
  \item{maxsize}{
    \emph{Maximum} physical size of the circles/squares plotted
    when \code{x} is a marked point pattern with 
    numerical marks.
    Incompatible with \code{meansize} and \code{markscale}.
  }
  \item{meansize}{
    \emph{Average} physical size of the circles/squares plotted
    when \code{x} is a marked point pattern with 
    numerical marks.
    Incompatible with \code{maxsize} and \code{markscale}.
  }
  \item{markscale}{
    physical scale factor determining the sizes of the
    circles/squares plotted when \code{x} is a marked point pattern with 
    numerical marks. Mark value will be multiplied by \code{markscale}
    to determine physical size.
    Incompatible with \code{maxsize} and \code{meansize}.
  }
  \item{minsize}{
    \emph{Minimum} physical size of the circles/squares plotted
    when \code{x} is a marked point pattern with 
    numerical marks.
    Incompatible with \code{zerosize}.
  }
  \item{zerosize}{
    Physical size of the circle/square representing a mark value of zero,
    when \code{x} is a marked point pattern with numerical marks.
    Incompatible with \code{minsize}.
    Defaults to zero.
  }
  \item{marktransform}{
    Experimental. A function that should be applied to the mark values
    before the symbol mapping is applied.
  }
}
\details{
  This algorithm determines a symbol map that can be used to 
  represent the points of \code{x} graphically.
  It serves as the default symbol map for the plot method 
  \code{\link{plot.ppp}}.

  Users can modify the behaviour of \code{\link{plot.ppp}} by
  saving the symbol map produced by \code{default.symbolmap},
  modifying the symbol map using \code{\link{update.symbolmap}} or other
  tools, and passing the modified symbol map to \code{\link{plot.ppp}}
  as the argument \code{symap}.
  
  The default representation depends on the marks of the points, as follows.
  \describe{
    \item{unmarked point pattern:}{
      If the point pattern does not have marks,
      then every point will be represented by the same plot symbol.
    }
    \item{multitype point pattern:}{
      If \code{marks(x)} is a factor, then 
      each level of the factor is 
      represented by a different plot character.
    }
    \item{continuous marks:}{
      If \code{marks(x)} is a numeric vector,
      each point is represented by a circle 
      with \emph{diameter} proportional to the mark
      (if the value is positive) or a square with \emph{side length}
      proportional to the absolute value of the mark
      (if the value is negative).
    }
    \item{other kinds of marks:}{
      If \code{marks(x)} is neither numeric nor a factor,
      then each possible mark will be represented by a
      different plotting character. The default is to
      represent the \eqn{i}th smallest mark value by
      \code{points(..., pch=i)}.
    }
  }
  The following arguments can be used to modify how the points are plotted:
  \itemize{
    \item
    If \code{fixsize=TRUE}, or if the graphics parameter \code{size} is
    given and is a single value, then numerical marks will be rendered
    as symbols of the same physical size
    \item
    The argument \code{chars} determines the plotting character
    or characters used to display the points (in all cases except
    for the case of continuous marks). For an unmarked point pattern,
    this should be a single integer or character determining a
    plotting character (see \code{par("pch")}).
    For a multitype point pattern, \code{chars} should be a vector
    of integers or characters, of the same length
    as \code{levels(marks(x))}, and then the \eqn{i}th level or type 
    will be plotted using character \code{chars[i]}.
    \item
    If \code{chars} is absent, but there is an extra argument
    \code{pch}, then this will determine the plotting character for
    all points.
    \item
    The argument \code{cols} determines the colour or colours used to
    display the points. For an unmarked point pattern,
    \code{cols} should be a character string
    determining a colour. For a multitype point pattern, \code{cols}
    should be a character vector, of the same length
    as \code{levels(marks(x))}: that is, there is one colour for each
    possible mark value. The \eqn{i}th level or type will
    be plotted using colour \code{cols[i]}. For a point pattern with
    continuous marks, \code{cols} can be either a character string
    or a character vector specifying colour values: the range of mark
    values will be mapped to the specified colours.
    Alternatively, for any kind of data, \code{cols} can be a colour map
    (object of class \code{"colourmap"}) created by \code{\link{colourmap}}.
    \item
    If \code{cols} is absent, the colours used to plot the
    points may be determined by the extra arguments \code{fg}
    and \code{bg} for foreground (edge) and background (fill) colours. 
    (These parameters are not recommended for plotting
    multitype point patterns, due to quirks of the graphics system.)
    \item
    The default colour for the points is a semi-transparent grey,
    if this is supported by the plot device. This behaviour can be
    suppressed (so that the default colour is non-transparent)
    by setting \code{spatstat.options(transparent=FALSE)}.
    \item 
    The arguments \code{maxsize}, \code{meansize} and \code{markscale} 
    are incompatible with each other (and incompatible with
    \code{symap}).
    The arguments \code{minsize} and \code{zerosize} are incompatible
    with each other (and incompatible with \code{symap}).
    Together, these arguments control the physical size of the circles and
    squares which represent the marks in a point pattern with continuous
    marks. The size of a circle is defined as its \emph{diameter};
    the size of a square is its side length.
    If \code{markscale} is given, then a mark value of \code{m}
    is plotted as a circle of diameter \code{m * markscale + zerosize}
    (if \code{m} is positive) or a square of side
    \code{abs(m) * markscale + zerosize}
    (if \code{m} is negative). If \code{maxsize} is given, then the
    largest mark in absolute value, \code{mmax=max(abs(marks(x)))},
    will be scaled to have physical size \code{maxsize}.
    If \code{meansize} is given, then the
    average absolute mark value, \code{mmean=mean(abs(marks(x)))},
    will be scaled to have physical size \code{meansize}.
    If \code{minsize} is given, then the
    minimum mark value, \code{mmean=mean(abs(marks(x)))},
    will be scaled to have physical size \code{minsize}.
    \item
    The user can set the default values of these plotting parameters
    using \code{\link{spatstat.options}("par.points")}.
  }

  Additionally the user can specify any of the graphics parameters
  recognised by \code{\link{symbolmap}}, including
  \code{shape}, \code{size}, \code{pch}, \code{cex},
  \code{cols}, \code{col}, \code{fg}, \code{bg}, 
  \code{lwd}, \code{lty}, \code{etch},
  \code{direction},\code{headlength},\code{headangle},\code{arrowtype}.

}
\value{
  A symbol map (object of class \code{"symbolmap"})
  or a list of symbol maps, one for each column of marks.
}
\author{
  \adrian
}
\seealso{
  \code{\link{plot.ppp}}

  \code{\link{default.symbolmap}}

  \code{\link{symbolmap}}
}
\examples{
  default.symbolmap(longleaf)
  default.symbolmap(lansing)
}
\keyword{spatial}
\keyword{hplot}
