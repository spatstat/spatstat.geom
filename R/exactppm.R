#'  exactppm.R
#'
#'  An internal device to represent Poisson point process models
#'  for which the MLE is computable exactly.
#'       - uniform intensity
#'       - intensity proportional to baseline
#'  These models are used mainly as a mathematical device
#'  in nonparametric methods such as 'rhohat' to represent the null/reference model,
#'  so that the code for nonparametric methods does not depend on 'ppm'
#'
#'  $Revision: 1.4 $ $Date: 2022/05/20 05:15:48 $

exactppm <- function(X, baseline=NULL, ..., subset=NULL, eps=NULL, dimyx=NULL) {
  stopifnot(inherits(X, c("ppp", "quad")))
  if(is.quad(X)) X <- X$data
  marx <- marks(X) # may be null
  lev <- levels(marx) # may be null
  if(is.null(subset)) {
    Xfit <- X
  } else {
    verifyclass(subset, "owin")
    Xfit <- X[subset]
  }
  if(is.null(baseline)) {
    #' stationary Poisson process
    beta <- intensity(Xfit)
  } else {
    #' Poisson process with intensity proportional to baseline
    if(is.im(baseline)) {
      denom <- integral(baseline, domain=subset)
    } else if(is.imlist(baseline) &&
              is.multitype(X) &&
              length(baseline) == length(lev)) {
      denom <- sapply(baseline, integral, domain=subset)
    } else if(is.function(baseline)) {
      if(!is.multitype(X) || length(formals(baseline)) == 2) {
        ba <- as.im(baseline, W=Window(X), eps=eps, dimyx=dimyx)
        denom <- integral(ba, domain=subset)
      } else {
        ba <- lapply(lev,
                     function(z) { as.im(baseline, W=Window(X), z, eps=eps, dimyx=dimyx)})
        denom <- sapply(ba, integral, domain=subset)
      }
    } else if(identical(baseline, "x")) {
      ba <- as.im(function(x,y){x}, W=Window(X), eps=eps, dimyx=dimyx)
      denom <- integral(ba, domain=subset)
    } else if(identical(baseline, "y")) {
      ba <- as.im(function(x,y){y}, W=Window(X), eps=eps, dimyx=dimyx)
      denom <- integral(ba, domain=subset)
    } else if(is.numeric(baseline) &&
              (length(baseline) == 1 ||
               is.multitype(X) && length(baseline) == length(lev))) {
      denom <- baseline * area(Window(Xfit))
    } else stop("Format of 'baseline' is not understood")
    numer <- if(is.multitype(Xfit)) as.integer(table(marks(Xfit))) else npoints(Xfit)
    beta <- numer/denom
    if(length(beta) == length(lev))
      names(beta) <- lev
  }
  model <- list(X=X, baseline=baseline, subset=subset, beta=beta)
  class(model) <- c("exactppm", class(model))
  return(model)
}

print.exactppm <- function(x, ...) {
  with(x, {
    splat("Exactly-fitted point process model")
    if(!is.multitype(X)) {
      if(is.null(baseline)) splat("Homogeneous intensity", signif(beta, 4)) else
                            splat("Intensity proportional to baseline", 
                                  paren(paste("proportionality constant",
                                              signif(beta, 4))))
    } else {
      lab <- levels(marks(X))
      if(is.null(baseline)) {
        splat("Homogeneous intensities:")
        splat(paste(paste(lab, signif(beta, 4), sep=": "), collapse=", "))
      } else {
        splat("Intensities proportional to baseline")
        splat("Proportionality constants:")
        splat(paste(paste(lab, signif(beta, 4), sep=": "), collapse=", "))
      }
    }
  })
  return(invisible(NULL))
}

predict.exactppm <- function(object, locations=NULL, ..., eps=NULL, dimyx=NULL) {
  X        <- object$X
  beta     <- object$beta # numeric
  baseline <- object$baseline # covariate or NULL
  if(is.null(locations))
    locations <- Window(X)
  if(length(beta) > 1) {
    ## Intensities for different types
    ## Syntax of 'evaluateCovariate' requires a list in this case
    beta <- as.list(beta)
  }
  ## evaluate at desired locations
  Beta <- evaluateCovariate(beta, locations, eps=eps, dimyx=dimyx)
  if(is.null(baseline)) {
    Lambda <- Beta
  } else {
    Baseline <- evaluateCovariate(baseline, locations, eps=eps, dimyx=dimyx)
    if(is.im(Beta) || is.imlist(Beta) || is.im(Baseline) || is.imlist(Baseline)) {
      Lambda <- imagelistOp(Beta, Baseline, "*")
    } else {
      Lambda <- Beta * Baseline
    }
  }
  ## tidy
  if(is.imlist(Lambda)) {
    if(length(Lambda) == 1) { 
      Lambda <- Lambda[[1]]
    } else if(length(Lambda) == length(Beta)) {
      names(Lambda) <- names(beta)
    }
  }
  return(Lambda)
}

