#'   NAobjects.R
#'
#'   A missing/unavailable object of class 'foo'
#'   inherits class c('NAobject', 'foo')
#'
#'   Methods for class 'NAobject' capture dispatch of print, plot, summary
#'   so that we don't need to tinker with print.ppp, plot.ppp etc.
#' 
#'   $Revision: 1.7 $ $Date: 2025/09/08 06:04:07 $
#' 
#'   ------------------------------------------------------------
#'        Make an 'NA object' of any class
#'   ------------------------------------------------------------

NAobject <- function(cls) {
  check.1.string(cls)
  structure(list(input = NA_character_,
                 wkt = NA_character_),
            class=c("NAobject", cls))
}

#'   ------------------------------------------------------------
#'        Recognise any 'NA object'
#'   ------------------------------------------------------------

is.NAobject <- function(x) { inherits(x, "NAobject") }

#'   ------------------------------------------------------------
#'        Common idiom to extract class ignoring 'NAobject'
#'   ------------------------------------------------------------

classIgnoringNA <- function(x, first=FALSE) {
  a <- setdiff(class(x), "NAobject")
  if(first) a <- a[1L]
  return(a)
}

#'   ------------------------------------------------------------
#'       methods for class 'NAobject'
#'   ------------------------------------------------------------

plot.NAobject <- function(x, ...) {
  splat("NA object: nothing plotted")
  invisible(NULL)
}

print.NAobject <- function(x, ...) {
  oc <- classIgnoringNA(x, first=TRUE)
  splat(paste0("<NA ", oc, ">"))
  invisible(NULL)
}

summary.NAobject <- function(object, ...) {
  oc <- classIgnoringNA(object, first=TRUE)
  structure(list(print=paste("NA object of class", sQuote(oc))),
            class="summary.NAobject")
}

print.summary.NAobject <- function(x, ...) {
  splat(x$print)
  invisible(NULL)
}

as.data.frame.NAobject <- function(x, ...) { NAobject("data.frame") }

is.marked.NAobject <- function(X, ...) { NA }

is.multitype.NAobject <- function(X, ...) { NA }
