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
#'  $Revision: 1.2 $ $Date: 2022/05/19 06:56:34 $

exactppm <- function(X, baseline=NULL, ..., subset=NULL, eps=NULL, dimyx=NULL) {
  stopifnot(inherits(X, c("ppp", "quad")))
  if(is.quad(X)) X <- X$data
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
              length(baseline) == length(levels(marks(X)))) {
      denom <- sapply(baseline, integral, domain=subset)
    } else if(is.function(baseline)) {
      if(!is.multitype(X) || length(formals(baseline)) == 2) {
        ba <- as.im(baseline, W=Window(X), eps=eps, dimyx=dimyx)
        denom <- integral(ba, domain=subset)
      } else {
        ba <- lapply(levels(marks(X)),
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
               is.multitype(X) && length(baseline) == length(levels(marks(X))))) {
      denom <- baseline * area(Window(Xfit))
    } else stop("Format of 'baseline' is not understood")
    numer <- if(is.multitype(Xfit)) as.integer(table(marks(Xfit))) else npoints(Xfit)
    beta <- numer/denom
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
  X <- object$X
  beta <- object$beta
  baseline <- object$baseline
  if(is.null(locations)) locations <- Window(X)
  locations <- as.mask(locations, eps=eps, dimyx=dimyx)
  betaimages <- lapply(beta, as.im, W=locations)
  if(is.null(baseline)) {
    lambdaimages <- betaimages
  } else if(is.im(baseline)) {
    lambdaimages <- solapply(betaimages, "*", e2=baseline)
  } else if(is.imlist(baseline)) {
    lambdaimages <- as.solist(mapply("*", e1=betaimages, e2=baseline, SIMPLIFY=FALSE))
  } else if(is.function(baseline)) {
    if(!is.multitype(X) || length(formals(baseline)) == 2) {
      baz <- as.im(baseline, W=locations, eps=eps, dimyx=dimyx)
      lambdaimages <- solapply(betaimages, "*", e2=baz)
    } else {
      baz <- lapply(levels(marks(X)),
                    function(z) { as.im(baseline, W=locations, z, eps=eps, dimyx=dimyx)})
      lambdaimages <- as.solist(mapply("*", e1=betaimages, e2=baz, SIMPLIFY=FALSE))
    }
  } else if(identical(baseline, "x")) {
    bxx <- as.im(function(x,y){x}, W=Window(X), eps=eps, dimyx=dimyx)
    lambdaimages <- solapply(betaimages, "*", e2=bxx)
  } else if(identical(baseline, "y")) {
    byy <- as.im(function(x,y){y}, W=Window(X), eps=eps, dimyx=dimyx)
    lambdaimages <- solapply(betaimages, "*", e2=byy)
  } else if(is.numeric(baseline)) {
    if(length(baseline) == 1) {
      lambdaimages <- solapply(betaimages, "*", e2=baseline)
    } else {
      lambdaimages <- as.solist(mapply("*", e1=betaimages, e2=baseline, SIMPLIFY=FALSE))
    }
  } else stop("Format of 'baseline' is not understood")
  return(if(length(lambdaimages) == 1) lambdaimages[[1]] else lambdaimages)
}
