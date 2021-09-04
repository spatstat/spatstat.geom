\name{invoke.metric}
\alias{invoke.metric}
\title{
  Perform Geometric Task using a Specified Metric
}
\description{
  Perform a desired geometrical operation
  using a specified distance metric.
}
\usage{
  invoke.metric(m, task, \dots, evaluate=TRUE)
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
  \item{evaluate}{
    Logical value specifying whether to actually perform the
    computation and return the result (\code{evaluate=TRUE}, the default)
    or to simply return the function which performs the computation
    (\code{evaluate=FALSE}).
  }
}
\details{
  A \sQuote{metric} is a measure of distance between points in space.
  An object of class \code{"metric"} represents such a metric,
  and supports many geometrical computations that involve
  the metric. See \code{\link{metric.object}}.

  The argument \code{task} should be the name of an existing function
  in the \pkg{spatstat} family representing a geometrical operation, such
  as computing pairwise distances, nearest-neighbour distances,
  the distance map, and so on. The code will determine whether
  this geometrical operation has a counterpart using the specified metric,
  that is defined and supported in the object \code{m}.
  If so, then this operation will be applied to the
  data specified in \code{\dots}, and the result will be returned.

  For example, the \pkg{spatstat} function \code{\link{nndist.ppp}} computes
  nearest-neighbour distances using the Euclidean distance metric.
  To calculate nearest-neighbour distances for a point pattern \code{X}
  using another metric \code{m},
  use \code{invoke.metric(m, "nndist.ppp", X)}.

  If \code{evaluate=FALSE}, the computation is not performed,
  and \code{invoke.metric} simply returns a function to
  perform the desired operation. 
}
\value{
  If \code{evaluate=TRUE} (the default), the result of the computation has the
  same format as the result of the computation using the existing
  function named \code{task}.
  
  If \code{evaluate=FALSE}, the result is a \code{function} in the \R language
  to perform the desired operation;
  or \code{NULL} if the operation is not supported by the metric.
}
\author{
  \adrian.
}
\seealso{
  \code{\link{convexmetric}}
}
\examples{
  ## nearest-neighbour distances using rectangular metric (L^1 metric)
  d <- convexmetric(square(c(-1,1)))
  y <- invoke.metric(d, "nndist.ppp", cells)
  f <- invoke.metric(d, "nndist.ppp", cells, evaluate=FALSE)
  y <- f(cells)
  invoke.metric(d, "orderPizza", evaluate=FALSE)
}
\keyword{spatial}
\keyword{math}
\concept{metric}
