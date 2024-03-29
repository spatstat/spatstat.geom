\name{progressreport}
\alias{progressreport}
\title{Print Progress Reports}
\description{
  Prints Progress Reports during a loop or iterative calculation.
}
\usage{
progressreport(i, n,
               every = min(100,max(1, ceiling(n/100))),
               tick = 1,
               nperline = NULL,
               charsperline = getOption("width"),
               style = spatstat.options("progress"),
               showtime = NULL,
               state=NULL,
               formula = (time ~ i),
               savehistory=FALSE)
}
\arguments{
  \item{i}{
    Integer. The current iteration number (from 1 to \code{n}).
  }
  \item{n}{
    Integer. The (maximum) number of iterations to be computed.
  }
  \item{every}{
    Optional integer. Iteration number will be printed 
    when \code{i} is a multiple of \code{every}.
  }
  \item{tick}{
    Optional integer. A tick mark or dot will be printed
    when \code{i} is a multiple of \code{tick}.
  }
  \item{nperline}{
    Optional integer. Number of iterations per line of output.
  }
  \item{charsperline}{
    Optional integer. The number of characters in a line of output.
  }
  \item{style}{
    Character string determining the style of display.
    Options are \code{"tty"} (the default), \code{"tk"} and \code{"txtbar"}.
    See Details.
  }
  \item{showtime}{
    Optional. Logical value indicating whether to print the estimated
    time remaining. Applies only when \code{style="tty"}.
  }
  \item{state}{
    Optional. A list containing the internal data.
  }
  \item{formula}{
    Optional. A model formula expressing the expected relationship between
    the iteration number \code{i} and the clock time \code{time}. Used for
    predicting the time remaining. 
  }
  \item{savehistory}{
    Optional. Logical value indicating whether to save the 
    elapsed times at which \code{progressreport} was called.
  }
}
\details{
  This is a convenient function for reporting progress
  during an iterative sequence of calculations
  or a suite of simulations.

  \itemize{
    \item 
    If \code{style="tk"} then \code{tcltk::tkProgressBar} is
    used to pop-up a new graphics window showing a progress bar.
    This requires the package \pkg{tcltk}.
    As \code{i} increases from 1 to \code{n}, the bar will lengthen.
    The arguments \code{every, tick, nperline, showtime} are ignored.
    \item 
    If \code{style="txtbar"} then \code{\link[utils]{txtProgressBar}} is
    used to represent progress as a bar made of text characters in the
    \R interpreter window.
    As \code{i} increases from 1 to \code{n}, the bar will lengthen.
    The arguments \code{every, tick, nperline, showtime} are ignored.
    \item 
    If \code{style="tty"} (the default),
    then progress reports are printed to the
    console. This only seems to work well under Linux.
    As \code{i} increases from 1 to \code{n}, 
    the output will be a sequence of dots (one dot for every \code{tick}
    iterations), iteration numbers (printed when iteration number is
    a multiple of \code{every} or is less than 4),
    and optionally the estimated time remaining
    and the estimated completion time.

    The estimated time remaining will be printed only if
    \code{style="tty"}, and the argument \code{state} is given,
    and either \code{showtime=TRUE}, or \code{showtime=NULL} and the
    iterations are slow (defined as: the estimated time remaining
    is longer than 3 minutes, or the average time per iteration is
    longer than 20 seconds).

    The estimated completion time will be printed only if the
    estimated time remaining is printed and the remaining time is
    longer than 10 minutes.

    By default, the estimated time remaining is calculated 
    by assuming that each iteration takes the same amount of time,
    and extrapolating.
    Alternatively, if the argument \code{formula} is given,
    then it should be a model formula, stating the expected relationship
    between the iteration number \code{i} and the clock time
    \code{time}. This model will be fitted to the history of clock times
    recorded so far, and used to predict the time remaining.
    (The default formula states that clock time is a linear function of the
    iteration number, which is equivalent to assuming that each
    iteration takes the same amount of time.)
  }
  
  It is optional, but strongly advisable, to use the argument \code{state}
  to store and update the internal data for the progress reports
  (such as the cumulative time taken for computation)
  as shown in the last example below.
  This avoids conflicts with other programs that might be
  calling \code{progressreport} at the same time.
}
\value{
  If \code{state} was \code{NULL}, the result is \code{NULL}.
  Otherwise the result is the updated value of \code{state}.
}
\author{
  \spatstatAuthors.
}
\examples{
  for(i in 1:40) {
     #
     # code that does something...
     # 
     progressreport(i, 40)
  }

  # saving internal state: *recommended*
  sta <- list()
  for(i in 1:20) {
     # some code ...
     sta <- progressreport(i, 20, state=sta)
  }

  #' use text progress bar
  sta <- list()
  for(i in 1:10) {
     # some code ...
     sta <- progressreport(i, 10, state=sta, style="txtbar")
  }

}
\keyword{print}
