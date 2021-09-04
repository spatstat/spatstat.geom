#'
#'   distancemetrics.R
#'
#'   Metrics on the spatial domain
#' 
#'   $Revision: 1.6 $ $Date: 2021/09/04 10:01:54 $

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

invoke.metric <- function(m, task, ..., evaluate=TRUE) {
  verifyclass(m, "metric")
  check.1.string(task)
  j <- match(task, names(m$tasks))
  if(is.na(j)) return(NULL)
  fname <- m$tasks[[j]]
  f <- m$functions[[fname]]
  if(!evaluate)
    return(f)
  if(is.null(f))
    stop(paste("This metric does not support", sQuote(task)), call.=FALSE)
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

