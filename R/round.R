#
#   round.R
#
#   discretisation of coordinates
#
#   $Revision: 1.7 $  $Date: 2024/06/12 06:23:53 $

round.ppp <- round.pp3 <- round.ppx <- function(x, digits=0, ...) {
  coords(x) <- round(as.matrix(coords(x)), digits=digits, ...)
  return(x)
}

rounding.ppp <- rounding.pp3 <- rounding.ppx <- function(x) {
  rounding(as.matrix(coords(x)))
}

