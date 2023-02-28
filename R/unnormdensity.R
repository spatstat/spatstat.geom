#
#  unnormdensity.R
#
#  $Revision: 1.13 $  $Date: 2023/02/28 02:37:33 $
#

unnormdensity <- function(x, ..., weights=NULL, defaults=list()) {
  if(any(!nzchar(names(list(...)))))
    stop("All arguments must be named (tag=value)")
  if(is.null(weights)) {
    ## all weights are 1 (not 1/n)
    out <- do.call.matched(density.default,
                           c(list(x=quote(x), ...), defaults))
    out$y <- length(x) * out$y
  } else if(length(weights) == 1) {
    ## all weights are equal
    out <- do.call.matched(density.default,
                           c(list(x=quote(x), ...), defaults))
    out$y <- weights[1] * length(x) * out$y
  } else if(length(weights) != length(x)) {
    stop("'x' and 'weights' have unequal length")
  } else if(all(weights == 0)) {
    ## result is zero
    out <- do.call.matched(density.default,
                           c(list(x=quote(x), ...), defaults))
    out$y <- 0 * out$y
  } else if(all(weights >= 0)) {
    # all masses are nonnegative, some are positive
    out <- do.call.matched(density.default,
                           c(list(x=quote(x),
                                  weights=quote(weights),
                                  subdensity=TRUE,
                                  ...),
                             defaults))
  } else if(all(weights <= 0)) {
    # all masses are nonpositive, some are negative
    w <- (- weights)
    out <- do.call.matched(density.default,
                           c(list(x=quote(x),
                                  weights=quote(w),
                                  subdensity=TRUE,
                                  ...),
                             defaults))
    out$y <- (- out$y)
  } else {
    # mixture of positive and negative masses
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
                                defaults))
      bw <- dabs$bw
    }
    ## compute densities for positive and negative masses separately
    outpos <- do.call.matched(density.default,
                              resolve.defaults(list(x=quote(x),
                                                    bw=bw,
                                                    adjust=1,
                                                    weights=quote(wpos),
                                                    subdensity=TRUE,
                                                    ...),
                                               defaults,
                                               .StripNull=TRUE))
    outneg <- do.call.matched(density.default,
                              resolve.defaults(list(x=quote(x),
                                                    bw=bw,
                                                    adjust=1,
                                                    weights=quote(wneg),
                                                    subdensity=TRUE,
                                                    ...),
                                               defaults,
                                               .StripNull=TRUE))
    ## combine
    out <- outpos
    out$y <- outpos$y - outneg$y
  }
  out$call <- match.call()
  return(out)
}

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
