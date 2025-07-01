#'   NAobjects.R
#'
#'   A missing/unavailable object of class 'foo'
#'   inherits class c('NAobject', 'foo')
#'
#'   Methods for class 'NAobject' capture dispatch of print, plot, summary
#'   so that we don't need to tinker with print.ppp, plot.ppp etc.
#' 
#'   $Revision: 1.3 $ $Date: 2025/06/30 12:30:36 $
#' 
#'   ------------------------------------------------------------
#'               specific NA objects 
#'   ------------------------------------------------------------

NA_ppp_ <- structure(list(input = NA_character_,
                          wkt = NA_character_),
                     class=c("NAobject", "ppp"))

NA_im_ <- structure(list(input = NA_character_,
                          wkt = NA_character_),
                     class=c("NAobject", "im"))

NA_owin_ <- structure(list(input = NA_character_,
                          wkt = NA_character_),
                     class=c("NAobject", "owin"))


#'   ------------------------------------------------------------
#'        Recognise any 'NA object'
#'   ------------------------------------------------------------

is.NAobject <- function(x) { inherits(x, "NAobject") }


#'   ------------------------------------------------------------
#'      methods for is.na (use only when the class is known)
#'   ------------------------------------------------------------

is.na.NAobject <- function(x) { TRUE }  # will generally be dispatched first

is.na.ppp <- function(x) { identical(x, NA_ppp_) } 

is.na.im <- function(x) { identical(x, NA_im_) }

is.na.owin <- function(x) { identical(x, NA_owin_) }


#'   ------------------------------------------------------------
#'       methods for class 'NAobject'
#'   ------------------------------------------------------------

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

