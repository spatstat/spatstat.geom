#'
#'   breakptgeom.R
#'
#'   Functions for creating a 'breakpts' object
#'   that depend on geometry of window, etc.
#'
#'   This code was excised from 'breakpts.R'
#'
#'   handle.r.b.args      Determine breakpoints for use in summary functions
#'                        such as Kest, Gest, Fest which recognise
#'                        arguments 'r' and 'breaks' and for which the
#'                        defaults depend on window geometry.
#'
#'   check.finespacing    Verify that breakpoint spacing is sufficiently fine
#'                        to ensure validity of discrete approximation to
#'                        product integral etc.
#' 
#'   $Revision: 1.2 $ $Date: 2023/11/05 00:58:19 $

handle.r.b.args <- function(r=NULL, breaks=NULL, window, pixeps=NULL,
                            rmaxdefault=NULL) {
  if(!is.null(r) && !is.null(breaks))
    stop(paste("Do not specify both",
               sQuote("r"), "and", sQuote("breaks")))
  if(!is.null(breaks)) {
    breaks <- as.breakpts(breaks)
  } else if(!is.null(r)) {
    breaks <- breakpts.from.r(r)
  } else {
    #' determine rmax
    #' ignore infinite or NA values of rmaxdefault
    if(!isTRUE(is.finite(rmaxdefault)))
      rmaxdefault <- NULL
    rmax <- rmaxdefault %orifnull% diameter(Frame(window))
    #' determine spacing
    if(is.null(pixeps)) {
      pixeps <-
        if(is.mask(window)) min(window$xstep, window$ystep) else rmax/128
    }
    rstep <- pixeps/4
    breaks <- make.even.breaks(rmax, bstep=rstep)
  }
  return(breaks)
}

check.finespacing <- function(r, eps=NULL, win=NULL,
                              rmaxdefault = max(r), 
                              context="",
                              action=c("fatal", "warn", "silent"),
                              rname) {
  if(missing(rname)) rname <- short.deparse(substitute(r))
  action <- match.arg(action)
  if(is.null(eps)) {
    b <- handle.r.b.args(window=win, rmaxdefault=rmaxdefault)
    eps <- b$step
  }
  dr <- max(diff(r))
  if(dr > eps * 1.01) {
    whinge <- paste(context, "the successive", rname,
                    "values must be finely spaced:",
                    "given spacing =",
                    paste0(signif(dr, 5), ";"),
                    "required spacing <= ",
                    signif(eps, 3))
    switch(action,
           fatal = stop(whinge, call.=FALSE),
           warn = warning(whinge, call.=FALSE),
           silent = {})
    return(FALSE)
  }
  return(TRUE)
}

