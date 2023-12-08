#'  pointweights.R
#'
#'  get a valid vector of weights for a point pattern
#'
#'  Argument 'weights' is usually passed from a user-level function
#'  It may be:
#'        a numeric vector
#'        a numeric matrix or data frame (if dfok=TRUE)
#'        a single number
#'        a function(x,y)
#'        a pixel image
#'        an expression involving the coordinates and marks
#' 
#'  $Revision: 1.4 $ $Date: 2023/12/08 06:51:20 $

pointweights <- function(X, ..., weights=NULL, parent=NULL, dfok=FALSE) {
  if(is.null(weights)) return(NULL)
  nX <- npoints(X)
  ## evaluate weights
  if(is.numeric(weights) && is.vector(as.numeric(weights))) {
    if(length(weights) == 1) weights <- rep(weights, nX)
  } else if(is.im(weights)) {
    weights <- safelookup(weights, X) # includes warning if NA
  } else if(is.function(weights)) {
      weights <- weights(X$x, X$y)
  } else if(is.expression(weights)) {
    #' evaluate expression in data frame of coordinates and marks
    df <- as.data.frame(X)
    weights <- try(eval(weights, envir=df, enclos=parent))
    if(inherits(weights, "try-error"))
      stop("Unable to evaluate expression for weights", call.=FALSE)
  } else if(dfok && inherits(weights, c("matrix", "data.frame"))) {
    weights <- as.matrix(weights)
  } else
    stop(paste0("Argument 'weights' should be ",
                "a numeric vector, ",
                if(dfok) "matrix or data frame, " else NULL,
                "a function, an image, or an expression"),
         call.=FALSE)

  ## trivial weights are returned as NULL
  if(length(weights) == 0) return(NULL)

  ## validate
  if(dfok && !is.null(dim(weights))) {
    check.nmatrix(weights, nX, squarematrix=FALSE, mname="weights")
  } else {
    check.nvector(weights, nX, vname="weights")
  }
  return(weights)
}

