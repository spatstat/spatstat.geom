\name{invoke.metric}
\alias{invoke.metric}
\alias{metricOperation}
\title{
  Perform Geometric Task using a Specified Metric
}
\description{
  Perform a desired geometrical operation
  using a specified distance metric.
}
\usage{
  invoke.metric(m, task, \dots)

  metricOperation(m, task)
}
\arguments{
  \item{m}{
    Metric (object of class \code{"metric"})
  }
  \item{task}{
    Character string specifying the task.
    The name of a function that performs
    the desired operation for the Euclidean metric. 
  }
  \item{\dots}{
    Input to the function that performs the geometrical operation
    (matching the arguments of \code{task}).
  }
}
\details{
  A \sQuote{metric} is a measure of distance between points in space.
  An object of class \code{"metric"} represents such a metric,
  and supports many geometrical computations that involve
  the metric. See \code{\link{metric.object}}.

  The argument \code{task} should be the name of an existing function
  in the \pkg{spatstat} family representing a geometrical operation, such
  as measuring pairwise distances, nearest-neighbour distances,
  calculating a distance map, and so on. The code will determine whether
  this geometrical operation has a counterpart for the specified metric,
  that is defined and supported in the object \code{m}.
  If so, then  \code{invoke.metric} will apply this operation to the
  data specified in \code{\dots}, while \code{metricOperation} simply
  returns a function that performs this operation.

  For example, the \pkg{spatstat} function \code{\link{nndist.ppp}} computes
  nearest-neighbour distances using the Euclidean metric.
  To calculate nearest-neighbour
  distances using another metric \code{m}, first call
  \code{f <- metricOperation(m, "nndist.ppp")} to obtain the function
  that calculates the nearest-neighbour distances, and then apply
  \code{f} to your point pattern data.

  If the desired operation is generic,
  then \code{task} should be the name
  of the particular \emph{method} that is required. 
}
\value{
  \code{metricOperation} returns a \code{function} in the \R language
  to perform the operation; or \code{NULL} if the operation is not supported.

  \code{invoke.metric} returns the result of the geometrical operation
  applied to \code{\dots}. This has the same format as the result
  of the function named \code{task}.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{convexmetric}}
}
\examples{
  ## nearest-neighbour distances using rectangular metric (L^1 metric)
  d <- convexmetric(square(1))
  y <- invoke.metric(d, "nndist.ppp", cells)
  metricOperation(d, "orderPizza")
}
\keyword{spatial}
\keyword{math}
\concept{metric}
