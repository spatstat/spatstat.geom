#
#  unnormdensity.R
#
#  $Revision: 1.9 $  $Date: 2022/03/24 01:24:57 $
#

unnormdensity <- function(x, ..., weights=NULL, defaults=list()) {
  if(any(!nzchar(names(list(...)))))
    stop("All arguments must be named (tag=value)")
  if(is.null(weights)) {
    out <- do.call.matched(density.default, c(list(x=x), list(...), defaults))
    out$y <- length(x) * out$y
  } else if(length(weights) == 1) {
    ## all weights are equal
    out <- do.call.matched(density.default, c(list(x=x), list(...), defaults))
    out$y <- weights[1] * out$y
  } else if(length(weights) != length(x)) {
    stop("'x' and 'weights' have unequal length")
  } else if(all(weights == 0)) {
    ## result is zero
    out <- do.call.matched(density.default, c(list(x=x), list(...), defaults))
    out$y <- 0 * out$y
  } else if(all(weights >= 0)) {
    # all masses are nonnegative, some are positive
    w <- weights
    totmass <- sum(w)
    out <- do.call.matched(density.default,
                           c(list(x=x),
                             list(...),
                             list(weights=w/totmass),
                             defaults))
    out$y <- out$y * totmass
  } else if(all(weights <= 0)) {
    # all masses are nonpositive, some are negative
    w <- (- weights)
    totmass <- sum(w)
    out <- do.call.matched(density.default,
                           c(list(x=x),
                             list(...),
                             list(weights=w/totmass),
                             defaults))
    out$y <- out$y * (- totmass)
  } else {
    # mixture of positive and negative masses
    w <- weights
    wabs <- abs(w)
    wpos <- pmax.int(0, w)
    wneg <- - pmin.int(0, w)
    ## determine bandwidth using absolute masses
    dabs <- do.call.matched(density.default,
                            c(list(x=x),
                              list(...),
                              list(weights=wabs/sum(wabs)),
                              defaults))
    bw <- dabs$bw
    ## compute densities for positive and negative masses separately
    sumwpos <- sum(wpos)
    sumwneg <- sum(wneg)
    outpos <- do.call.matched(density.default,
                              resolve.defaults(list(x=quote(x)),
                                               list(bw=bw, adjust=1),
                                               list(weights=wpos/sumwpos),
                                               list(...),
                                               defaults,
                                               .StripNull=TRUE))
    outneg <- do.call.matched(density.default,
                              resolve.defaults(list(x=quote(x)),
                                               list(bw=bw, adjust=1),
                                               list(weights=wneg/sumwneg),
                                               list(...),
                                               defaults,
                                               .StripNull=TRUE))
    ## combine
    out <- outpos
    out$y <- sumwpos * outpos$y - sumwneg * outneg$y
  }
  out$call <- match.call()
  return(out)
}

