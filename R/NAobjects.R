#'   NAobjects.R
#'
#'   A missing/unavailable object of class 'foo'
#'   inherits class c('NAobject', 'foo')
#'
#'   Methods for class 'NAobject' capture dispatch of print, plot, summary
#' 
#'   $Revision: 1.3 $ $Date: 2025/06/30 12:30:36 $
#' 

NA_ppp_ <- structure(list(input = NA_character_,
                          wkt = NA_character_),
                     class=c("NAobject", "ppp"))

NA_im_ <- structure(list(input = NA_character_,
                          wkt = NA_character_),
                     class=c("NAobject", "im"))

is.na.NAobject <- function(x) { TRUE }

is.na.ppp <- function(x) { identical(x, NA_ppp_) }  # could probably be set to return FALSE

is.na.im <- function(x) { identical(x, NA_im_) }


plot.NAobject <- function(x, ...) {
  splat("NA object: nothing plotted")
  invisible(NULL)
}

print.NAobject <- function(x, ...) {
  splat("<NA", paste0(class(x)[-1], ">"))
  invisible(NULL)
}

summary.NAobject <- function(object, ...) {
  oc <- class(object)[-1]
  structure(list(print=paste("NA object of class", sQuote(oc))),
            class="summary.NAobject")
}

print.summary.NAobject <- function(x, ...) {
  splat(x$print)
  invisible(NULL)
}

