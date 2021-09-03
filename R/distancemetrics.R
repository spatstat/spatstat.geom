#'
#'   distancemetrics.R
#'
#'   Metrics on the spatial domain
#' 
#'   $Revision$ $Date$

#'  An object of class 'metric' is a metric
#'  containing
#'     $functions: named list of internal functions
#'     $tasks: named vector of possible tasks, mapped to names of $functions

print.metric <- function(x, ...) { x$print() }

invoke.metric <- function(m, task, ...) {
  verifyclass(m, "metric")
  stopifnot(is.character(task))
  stopifnot(length(task) == 1)
  j <- match(task, names(m$tasks))
  if(is.null(j))
    stop(paste("This metric does not support", sQuote(task)), call.=FALSE)
  fname <- m$tasks[[j]]
  f <- m$functions[[fname]]
  f(...)
}

#'  An object of class 'metricfun' is a function that creates a metric

print.metricfun <- function(x, ...) {
  anames <- names(formals(anames))
  splat(paste0("function", paren(paste(anames,collapse=", "))))
  if(!is.null(ex <- attr(x, "explain")))
    splat(ex)
  return(invisible(NULL))
}

