#'
#'   distancemetrics.R
#'
#'   Metrics on the spatial domain
#' 
#'   $Revision: 1.5 $ $Date: 2021/09/04 04:29:02 $

#'  An object of class 'metric' is a metric
#'  containing
#'     $functions: named list of internal functions
#'     $tasks: named vector of possible tasks, mapped to names of $functions

print.metric <- function(x, ...) { x$functions$print() }

summary.metric <- function(object, ...) {
  print(object, ...)
  splat("Supported operations:", commasep(sQuote(names(object$tasks))))
  invisible(NULL)
}

invoke.metric <- function(m, task, ...) {
  f <- metricOperation(m, task)
  if(is.null(f))
    stop(paste("This metric does not support", sQuote(task)), call.=FALSE)
  f(...)
}

metricOperation <- function(m, task) {
  verifyclass(m, "metric")
  check.1.string(task)
  j <- match(task, names(m$tasks))
  if(is.na(j)) return(NULL)
  fname <- m$tasks[[j]]
  f <- m$functions[[fname]]
  return(f)
}

#'  An object of class 'metricfun' is a function that creates a metric

print.metricfun <- function(x, ...) {
  anames <- names(formals(anames))
  splat(paste0("function", paren(paste(anames,collapse=", "))))
  if(!is.null(ex <- attr(x, "explain")))
    splat(ex)
  return(invisible(NULL))
}

