#'
#'	breakpts.R
#'
#'	A simple class definition for the specification
#'      of histogram breakpoints for nonnegative numbers (such as distances)
#'      in the special form we need them.
#'
#'
#'	$Revision: 1.31 $	$Date: 2023/11/05 01:02:04 $
#'
#'      The breakpoints must
#'           (a) span the range of the data
#'           (b) be given in increasing order
#'           (c) satisfy breaks[2] = 0.
#'
#'	The function make.even.breaks() will create suitable breakpoints.
#'
#'      Condition (c) means that the first histogram cell has
#'      *right* endpoint equal to 0.
#'
#'      Since the numerical data are nonnegative, the effect of (c) is
#'      that the first histogram cell counts the number of values which are
#'      exactly equal to 0. Hence F(0), the probability P{X = 0},
#'      is estimated without a discretisation bias.
#'
#'	We assume the histograms have followed the default counting rule
#'	in hist.default(), which is such that the k-th entry of the histogram
#'	counts the number of data values in 
#'		I_k = ( breaks[k],breaks[k+1] ]	for k > 1
#'		I_1 = [ breaks[1],breaks[2]   ]
#'
#'	The implementations of estimators of distance distributions
#'      in the spatstat package return vectors of length = length(breaks)-1
#'      with value[k] = estimate of F(breaks[k+1]),
#'      i.e. value[k] is an estimate of the c.d.f. at the RIGHT endpoint
#'      of the kth histogram cell.
#'
#'      An object of class 'breakpts' contains:
#'
#'             $val     the actual breakpoints
#'             $max     the maximum value (= last breakpoint)
#'             $ncells  total number of histogram cells
#'             $r       right endpoints, r = val[-1]
#'             $even    logical = TRUE if cells known to be evenly spaced
#'             $npos    number of histogram cells on the positive halfline
#'                       = length(val) - 2,
#'                      or NULL if cells not evenly spaced
#'             $step    histogram cell width
#'                      or NULL if cells not evenly spaced
#'      
#' --------------------------------------------------------------------
#' 

breakpts <- function(val, maxi, even=FALSE, npos=NULL, step=NULL) {
  out <- list(val=as.numeric(val),
              max=as.numeric(maxi),
              ncells=length(val)-1L, r = val[-1L],
              even=isTRUE(even),
              npos=npos, step=step)
  class(out) <- "breakpts"
  out
}

make.even.breaks <- function(bmax, npos, bstep) {
  bmax <- as.numeric(bmax)
  if(bmax <= 0)
    stop("bmax must be positive")
  if(missing(bstep) && missing(npos))
    stop(paste("Must specify either", sQuote("bstep"),
               "or", sQuote("npos")))
  if(!missing(npos)) {
    npos <- as.integer(npos)
    bstep <- bmax/npos
    val <- seq(from=0, to=bmax, length.out=npos+1L)
    val <- c(-bstep,val)
    right <- bmax
  } else {
    bstep <- as.numeric(bstep)
    npos <- ceiling(bmax/bstep)
    right <- bstep * npos
    val <- seq(from=0, to=right, length.out=npos+1L)
    val <- c(-bstep,val)
  }
  breakpts(val, right, TRUE, npos, bstep)
}

as.breakpts <- function(...) {

  XL <- list(...)

  if(length(XL) == 1L) {
    #' There is a single argument
    X <- XL[[1L]]

    if(inherits(X, "breakpts")) {
      ## X already in correct form
      return(X)
    }
  
    if(is.vector(X) && length(X) > 2) {
      ## it's a vector
      X <- as.numeric(X)
      if(X[2L] != 0)
        stop("breakpoints do not satisfy breaks[2] = 0")
      #'The following test for equal spacing is used in hist.default
      steps <- diff(X)
      if(diff(range(steps)) < 1e-07 * mean(steps)) {
        #'equally spaced
        return(breakpts(X, max(X), TRUE, length(X)-2, steps[1L]))
      } else {
        #'unknown spacing
        return(breakpts(X, max(X), FALSE))
      }
    }
    
  } else {
    #' There are multiple arguments.
    
    #' Exactly two arguments 
    if(length(XL) == 2)
      return(make.even.breaks(XL[[1L]], XL[[2L]]))

    #' Two arguments named 'max' and 'npos'
  
    if(!is.null(XL$max) && !is.null(XL$npos))
      return(make.even.breaks(XL$max, XL$npos))

    #' otherwise
    stop("Don't know how to convert these data to breakpoints")
  }
  #' never reached
}


check.hist.lengths <- function(hist, breaks) {
  #' internal check for consistency between histogram result and breakpoints
  stopifnot(inherits(breaks, "breakpts"))
  nh <- length(hist)
  nb <- breaks$ncells
  if(nh != nb)
    stop(paste("Length of histogram =", nh,
               "not equal to number of histogram cells =", nb))
}

breakpts.from.r <- function(r) {
  if(!is.numeric(r) && !is.vector(r))
    stop("r must be a numeric vector")
  r <- as.numeric(r)
  if(length(r) < 2)
    stop(paste("r has length", length(r), "- must be at least 2"))
  if(r[1L] != 0)
    stop("First r value must be 0")
  if(any(diff(r) <= 0))
    stop("successive values of r must be increasing")
  dr <- r[2L] - r[1L]
  b <- c(-dr, r)
  return(as.breakpts(b))
}

