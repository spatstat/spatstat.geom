#'
#' evaluatecovariates.R
#'
#' Evaluate covariates at specified locations
#'
#'   $Revision: 1.8 $  $Date: 2022/05/20 06:25:05 $
#'

evaluateCovariate <- function(covariate, locations, ...) {
  if(is.owin(locations)) {
    evaluateCovariateAtPixels(covariate, locations, ...)
  } else {
    evaluateCovariateAtPoints(covariate, locations, ...)
  }
}

evaluateCovariateAtPoints <- function(covariate, locations, ...,
                                      allow.column=TRUE) {
  AvoidNames <- c("eps", "dimyx", "types")
  stopifnot(is.ppp(locations))
  n <- npoints(locations)
  marx <- marks(locations) # may be null
  lev <- levels(marx) # may be null
  if(is.im(covariate)) {
    ## single pixel image
    values <- safelookup(covariate, locations)
  } else if(is.imlist(covariate)) {
    ## list of images for each type of point
    if(length(covariate) != length(lev)) 
      stop(paste("Number of covariate images", paren(length(covariate)),
                 "does not match number of possible types in data",
                 paren(length(lev))),
           call.=FALSE)
    values <- vector(mode="list", length=n)
    mm <- as.integer(marx)
    for(k in which(as.integer(table(mm)) > 0)) {
      relevant <- which(mm == k)
      values[relevant] <- safelookup(covariate[[k]], locations[relevant])
    }
    values <- unlist(values)
  } else if(is.function(covariate)) {
    ## function(x,y) or function(x,y,m)
    if(length(formals(covariate)) <= 2L ||
       !any(c("m", "marks") %in% names(formals(covariate)))) {
      ## function does not use mark value
      values <- do.call.without(covariate,
                                locations$x, locations$y, ...,
                                avoid=AvoidNames)
    } else {
      ## function expects the mark values
      values <- do.call.without(covariate,
                                locations$x, locations$y, marx, ...,
                                avoid=AvoidNames)
    }
  } else if(is.list(covariate) && all(sapply(covariate, is.function))) {
    ## list of functions for each type of point
    if(length(covariate) != length(lev))
      stop(paste("Number of covariate functions", paren(length(covariate)),
                 "does not match number of possible types in data",
                 paren(length(lev))),
           call.=FALSE)
    values <- vector(mode="list", length=n)
    xx <- locations$x
    yy <- locations$y
    mm <- as.integer(marx)
    for(k in which(as.integer(table(mm)) > 0)) {
      relevant <- which(mm == k)
      values[relevant] <- do.call.without(covariate[[k]],
                                          xx[relevant], yy[relevant], ...,
                                          avoid=AvoidNames)
    }
    values <- unlist(values)
  } else if(is.numeric(covariate) || is.factor(covariate)) {
    ## numeric/categorical value or vector
    if(length(covariate) == 1L) {
      values <- rep.int(covariate, n)
    } else if(allow.column && length(covariate) == n) {
      ## NOTE
      values <- covariate 
    } else stop("Inappropriate length for covariate vector")
  } else if(is.list(covariate) && all(lengths(covariate) == 1L) &&
            (all(sapply(covariate, is.numeric)) || all(sapply(covariate, is.factor)))) {
    ## list of single values, assumed to correspond to types
    if(length(covariate) != length(lev)) 
      stop(paste("Length of list of covariate values", paren(length(covariate)),
                 "does not match number of possible types in data",
                 paren(length(lev))),
           call.=FALSE)
    values <- unlist(covariate[as.integer(marx)])
  } else 
    stop("Covariate should be an image, a function or a factor/numeric vector")
  return(values)
}

evaluateCovariateAtPixels <- function(covariate, locations, ...,
                                  types=NULL, eps=NULL, dimyx=NULL) {
  stopifnot(is.owin(locations))
  M <- as.mask(locations, eps=eps, dimyx=dimyx)
  if(is.im(covariate)) {
    value <- covariate[M, raster=M, drop=FALSE]
  } else if(is.imlist(covariate)) {
    value <- solapply(covariate, "[", i=M, raster=M, drop=FALSE)
  } else if(is.function(covariate)) {
    if(is.null(types) || length(formals(covariate)) <= 2L ||
       !any(c("m", "marks") %in% names(formals(covariate)))) {
      ## function (x,y) or function(x,y,...)  does not use mark value
      value <- as.im(covariate, W=M, ...)
    } else {
      ## function f(x,y,m) or function(x,y,m, ...) expects mark value
      value <- solapply(types, function(a) { as.im(covariate, W=M, a, ...) })
    }
  } else if(is.list(covariate) && all(sapply(covariate, is.function))) {
    ## list of function(x,y) or list of function(x,y,..)
    value <- solapply(covariate, as.im, W=M, ...)
  } else if(identical(covariate, "x")) {
    value <- as.im(function(x,y){x}, W=M)
  } else if(identical(covariate, "y")) {
    value <- as.im(function(x,y){y}, W=M)
  } else if(is.numeric(covariate) || is.factor(covariate)) {
    value <- solapply(covariate, as.im, W=M)
  } else if(is.list(covariate) && all(lengths(covariate) == 1L) &&
            all(sapply(covariate, is.numeric) | sapply(covariate, is.factor))) {
    ## list of single values, associated with types
    if(length(types) > 0 && length(covariate) != length(types))
      stop(paste("Length of list of covariate values", paren(length(covariate)),
                 "does not match number of possible types in data",
                 paren(length(types))),
           call.=FALSE)
    value <- solapply(covariate, as.im, W=M)
  } else stop("Format of covariate is not understood")
  if(!is.null(types)) {
    ## Ensure result is a solist of the right length
    if(is.im(value)) {
      value <- rep(list(value), length(types))
    } else if(length(value) != length(types)) 
      warning("Mismatch between number of covariates and number of types", call.=FALSE)
    names(value) <- types
  }
  return(value)
}


