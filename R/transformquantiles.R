## transformquantiles.R

## probability integral transformation
## aka histogram equalisation
## aka transformation to uniformity

## $Revision: 1.1 $ $Date: 2023/11/04 04:39:11 $

transformquantiles <- function(X, uniform=FALSE, reverse=FALSE, ...) {
  if(!uniform && !reverse) return(X)
  o <- order(X[])
  V <- X
  n <- length(o)
  if(uniform) V[][o] <- (seq_len(n) - 0.5)/n
  if(reverse) V[][o] <- V[][rev(o)]
  return(V)
}
