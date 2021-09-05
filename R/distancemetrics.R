#'
#'   distancemetrics.R
#'
#'   Metrics on the spatial domain
#' 
#'   $Revision: 1.9 $ $Date: 2021/09/05 03:05:47 $
#'
#'  An object of class 'metric' is essentially a named list of functions
#'  where the names specify the tasks. 
#'
#'  An object of class 'metricfun' is a function that creates a metric
#'
#'  See 'convexdist.R' for an example.

## ..............  metric ................................

print.metric <- function(x, ...) { x$print() }

summary.metric <- function(object, ...) {
  print(object, ...)
  splat("Supported operations:", commasep(sQuote(names(object))))
  invisible(NULL)
}

invoke.metric <- function(m, task, ..., evaluate=TRUE) {
  verifyclass(m, "metric")
  check.1.string(task)
  j <- match(task, names(m))
  f <- if(is.na(j)) NULL else m[[j]]
  if(!evaluate)
    return(f)
  if(is.null(f))
    stop(paste("This metric does not support", sQuote(task)), call.=FALSE)
  f(...)
}

## ..............  metricfun .............................

#'  An object of class 'metricfun' is a function that creates a metric

print.metricfun <- function(x, ...) {
  anames <- names(formals(x))
  splat(paste0("function", paren(paste(anames,collapse=", "))))
  if(!is.null(ex <- attr(x, "explain")))
    splat(ex)
  return(invisible(NULL))
}


## ......... Utilities to trap user errors ........................

## Utility for existing functions which do not support non-Euclidean metric)

warn.no.metric.support <- function(caller, ..., metric) {
  if(!missing(metric))
    warning(paste("Argument 'metric' is not implemented for",
                  paste0(sQuote(caller), " and was ignored")))
  invisible(NULL)
}

## Utility for use in metric counterparts of standard functions,
## when some arguments of standard function are unsupported by metric function
## (Issues a message only if the arguments have non-default values)

warn.unsupported.args <- function(unsup, ...) {
  given <- list(...)
  if(any(hit <- names(unsup) %in% names(given))) {
    values <- resolve.defaults(given, unsup)[names(unsup)]
    changed <- !mapply(identical, x=unsup, y=values)
    if(any(changed)) {
      n <- sum(changed)
      warning(paste(ngettext(n, "Argument", "Arguments"),
                    commasep(sQuote(names(unsup)[changed])),
                    ngettext(n, "is", "are"),
                    "not supported by this metric, and",
                    ngettext(n, "was", "were"), "ignored"),
              call.=FALSE)
    }
  }
  invisible(NULL)
}
