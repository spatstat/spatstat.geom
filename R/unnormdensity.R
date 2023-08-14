#
#  unnormdensity.R
#
#  $Revision: 1.19 $  $Date: 2023/08/14 06:27:35 $
#

unnormdensity <- local({

  unnormdensity <- function(x, ..., weights=NULL, defaults=list()) {
    if(any(!nzchar(names(list(...)))))
      stop("All arguments must be named (tag=value)")
    envir.here <- sys.frame(sys.nframe())
    ## suppress annoying warnings in density.default
    defaults <- resolve.defaults(defaults,
                                 list(warnWbw=FALSE))
    if(length(x) <= 1) {
      ## density.default does not handle this 
      out <- do.call(fewdatacase,
                     resolve.defaults(
                       list(x=x, ..., weights=weights),
                       defaults))
    } else if(is.null(weights)) {
      ## all weights are 1 (not 1/n)
      out <- do.call.matched(density.default,
                             c(list(x=quote(x), ...), defaults),
                             envir=envir.here)
      out$y <- length(x) * out$y
    } else if(length(weights) == 1) {
      ## all weights are equal
      out <- do.call.matched(density.default,
                             c(list(x=quote(x), ...), defaults),
                             envir=envir.here)
      out$y <- weights[1] * length(x) * out$y
    } else if(length(weights) != length(x)) {
      stop("'x' and 'weights' have unequal length")
    } else {
      weightrange <- range(weights)
      if(all(weightrange == 0)) {
        ## result is zero
        out <- do.call.matched(density.default,
                               c(list(x=quote(x), ...), defaults),
                               envir=envir.here)
        out$y <- 0 * out$y
      } else if(all(weightrange >= 0)) {
        ## all masses are nonnegative, some are positive
        out <- do.call.matched(density.default,
                               c(list(x=quote(x),
                                      weights=quote(weights),
                                      subdensity=TRUE,
                                      ...),
                                 defaults),
                               envir=envir.here)
      } else if(all(weightrange <= 0)) {
        ## all masses are nonpositive, some are negative
        w <- (- weights)
        out <- do.call.matched(density.default,
                               c(list(x=quote(x),
                                      weights=quote(w),
                                      subdensity=TRUE,
                                      ...),
                                 defaults),
                               envir=envir.here)
        out$y <- (- out$y)
      } else {
        ## mixture of positive and negative masses
        w <- weights
        wabs <- abs(w)
        wpos <- pmax.int(0, w)
        wneg <- - pmin.int(0, w)
        ## determine bandwidth value
        bw <- list(...)$bw  # could be NULL
        if(is.numeric(bw)) {
          ## bandwidth is given, as a numeric value
          ## adjust by factor 'adjust'
          adjust <- list(...)$adjust %orifnull% 1
          bw <- bw * adjust
        } else {
          ## compute bandwidth by applying a rule, using absolute masses
          dabs <- do.call.matched(density.default,
                                  c(list(x=quote(x),
                                         weights=quote(wabs),
                                         subdensity=TRUE,
                                         ...),
                                    defaults),
                                  envir=envir.here)
          bw <- dabs$bw
        }
        ## Bandwidth is now determined as a numeric value.
        ## Compute densities for positive and negative masses separately
        outpos <- do.call.matched(density.default,
                                  resolve.defaults(list(x=quote(x),
                                                        bw=bw,
                                                        adjust=1,
                                                        weights=quote(wpos),
                                                        subdensity=TRUE,
                                                        ...),
                                                   defaults,
                                                   .StripNull=TRUE),
                                  envir=envir.here)
        outneg <- do.call.matched(density.default,
                                  resolve.defaults(list(x=quote(x),
                                                        bw=bw,
                                                        adjust=1,
                                                        weights=quote(wneg),
                                                        subdensity=TRUE,
                                                        ...),
                                                   defaults,
                                                   .StripNull=TRUE),
                                  envir=envir.here)
        ## combine
        out <- outpos
        out$y <- outpos$y - outneg$y
      }
    }
    out$call <- match.call()
    return(out)
  }

  fewdatacase <- function(x, ...,
                          weights=NULL, kernel="gaussian",
                          bw=NULL,
                          from=NULL, to=NULL, n=512) {
    nx <- length(x)
    if(nx == 0) {
      needed <- list(from=from, to=to, n=n)
    } else if(nx == 1) {
      needed <- list(bw=bw, from=from, to=to, n=n)
    } else stop("Internal function 'fewdatacase' was invoked incorrectly")

    if(any(absent <- sapply(needed, is.null))) {
      nabsent <- sum(absent)
      stop(paste(ngettext(nabsent, "Argument", "Arguments"),
                 commasep(sQuote(names(needed)[absent])),
                 ngettext(nabsent, "is", "are"),
                 "required in density.default when length(x) <= 1"),
           call.=FALSE)
    }
    
    xx <- seq(from, to, length.out=n)

    if(nx == 0) {
      yy <- rep(0, n)
    } else {
      kernel <- match.kernel(kernel)
      if(!is.numeric(bw)) {
        ## Bandwidth selection rules require >= 2 data values
        ## Use fallback
        h <- (to-from)/3
        bw <- h / kernel.factor(kernel)
      }
      yy <- (weights %orifnull% 1) * dkernel(xx, kernel=kernel, mean=x, sd=bw)
    } 
    structure(list(x=xx, y=yy, bw=bw, n=n, call=match.call(),
                   data.name='x', has.na=FALSE),
              class="density")
  }

  unnormdensity
})


integral.density <- function(f, domain=NULL, weight=NULL, ...) {
  x <- f$x
  y <- f$y
  if(!is.null(domain)) {
    check.range(domain)
    retain <- inside.range(x, domain)
    x <- x[retain]
    y <- y[retain]
  }
  if(!is.null(weight)) {
    stopifnot(is.function(weight))
    y <- y * weight(x)
  }
  dx <- diff(x)
  ybar <- (y[-1] + y[-length(y)])/2
  sum(dx * ybar)
}
