#'  pointweights.R
#'
#'  get a valid vector of weights for a point pattern
#'
#'  (can also be used for evaluating a covariate at the data points)
#'
#'  Argument 'weights' is usually passed from a user-level function
#'  It may be:
#'        a numeric or logical vector
#'        a numeric or logical matrix or data frame (if dfok=TRUE)
#'        a single number
#'        a function(x,y)
#'        a pixel image
#'        an expression involving the coordinates and marks
#'        one of the characters "x" or "y"
#' 
#'  $Revision: 1.10 $ $Date: 2026/02/13 06:42:31 $

pointweights <- function(X, ..., weights=NULL, parent=NULL, dfok=FALSE,
                         weightsname="weights") {
  if(is.null(weights)) return(NULL)
  nX <- npoints(X)
  weightsblurb <- weightsname
  ## evaluate weights
  if(is.character(weights) && length(weights) == 1) {
    switch(weights,
           x = { return(X$x) },
           y = { return(X$y) },
           {
             ## Not a reserved symbol
             ## Treat as an expression and continue
             if(missing(weightsname)) weightsname <- weights
             weights <- str2expression(weights)
           })
  }
  if((is.numeric(weights) || is.logical(weights)) && NCOL(weights) == 1) {
    if(length(weights) == 1) weights <- rep(weights, nX)
    weights <- as.numeric(weights)
  } else if(is.im(weights)) {
    weights <- safelookup(weights, X) # includes warning if NA
    weightsblurb <- paste("Pixel values of", sQuote(weightsname))
  } else if(is.function(weights)) {
    weights <- try(weights(X$x, X$y))
    if(inherits(weights, "try-error"))
      stop(paste("Evaluating the function", sQuote(weightsname),
                 "returned an error"),
           call.=FALSE)
    weightsblurb <- paste("Values of the function", sQuote(weightsname))
  } else if(is.expression(weights)) {
    #' evaluate expression in data frame of coordinates and marks
    df <- as.data.frame(X)
    weights <- try(eval(weights, envir=df, enclos=parent))
    if(inherits(weights, "try-error"))
      stop(paste("Unable to evaluate expression for", weightsname), call.=FALSE)
    if(is.im(weights)) {
      weights <- safelookup(weights, X) # includes warning if NA
      weightsblurb <- paste("Pixel values of", sQuote(weightsname))
    } else {
      weightsblurb <- paste("Values of the expression", sQuote(weightsname))
    }
  } else if(NCOL(weights) > 1) {
    weights <- as.matrix(weights)
  } else
    stop(paste("Argument", sQuote(weightsname),
               "should be a numeric or logical",
               if(dfok) "vector, matrix or data frame, or" else "vector,",
                "a function, an image, or an expression"),
         call.=FALSE)

  ## trivial weights are returned as NULL
  if(length(weights) == 0) return(NULL)

  ## validate
  if(NCOL(weights) > 1) {
    if(!dfok) stop("A matrix or data frame of weights is not permitted",
                   call.=FALSE)
    weights <- apply(weights, 2, as.numeric)
    check.nmatrix(weights, nX, squarematrix=FALSE, mname=weightsblurb)
  } else {
    weights <- as.numeric(weights)
    check.nvector(weights, nX, vname=weightsblurb)
  }
  return(weights)
}

