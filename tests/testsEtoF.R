#'
#'   Header for all (concatenated) test files
#'
#'   Require spatstat.geom
#'   Obtain environment variable controlling tests.
#'
#'   $Revision: 1.5 $ $Date: 2020/04/30 05:31:37 $

require(spatstat.geom)
FULLTEST <- (nchar(Sys.getenv("SPATSTAT_TEST", unset="")) > 0)
ALWAYS   <- TRUE
cat(paste("--------- Executing",
          if(FULLTEST) "** ALL **" else "**RESTRICTED** subset of",
          "test code -----------\n"))
#  tests/emptymarks.R
#
# test cases where there are no (rows or columns of) marks
#
#  $Revision: 1.4 $ $Date: 2020/04/28 12:58:26 $

if(ALWAYS) {
local({
  n <- npoints(cells)
  df <- data.frame(x=1:n, y=factor(sample(letters, n, replace=TRUE)))
  nocolumns <- c(FALSE, FALSE)
  norows <- rep(FALSE, n)
  X <- cells
  marks(X) <- df
  marks(X) <- df[,1]
  marks(X) <- df[,nocolumns]
  Z <- Y <- X[integer(0)]
  marks(Y) <- df[norows,]
  stopifnot(is.marked(Y))
  marks(Z) <- df[norows,nocolumns]
  stopifnot(!is.marked(Z))
})
}
#
#    tests/factorbugs.R
#
# check for various bugs related to factor conversions
#
#    $Revision: 1.6 $  $Date: 2020/12/03 03:26:25 $
#

if(ALWAYS) {
local({
  ## make a factor image
  m <- factor(rep(letters[1:4], 4))
  Z <- im(m, xcol=1:4, yrow=1:4)
  ## make a point pattern
  set.seed(42)
  X <- runifrect(20, win=as.owin(Z))
  ## look up the image at the points of X
  ## (a) internal
  ans1 <- lookup.im(Z, X$x, X$y)
  stopifnot(is.factor(ans1))
  ## (b) user level
  ans2 <- Z[X]
  stopifnot(is.factor(ans2))
  ## (c) turn the image into a tessellation
  ##  and apply quadratcount
  V <- tess(image = Z)
  quadratcount(X, tess=V)
  ## (d) pad image
  Y <- padimage(Z, factor("b", levels=levels(Z)))
  stopifnot(Y$type == "factor")
  U <- padimage(Z, "b")
  stopifnot(U$type == "factor")
  ## (e) manipulate levels
  Zb <- relevel(Z, "b")
  Zv <- mergeLevels(Z, vowel="a", consonant=c("b","c","d"))
  P <- X %mark% Z[X]
  Pv <- mergeLevels(P, vowel="a", consonant=c("b","c","d"))
})
}
#
#  tests/func.R
#
#   $Revision: 1.9 $   $Date: 2022/10/23 00:48:40 $
#
#  Tests of 'funxy' infrastructure etc

if(FULLTEST) {
local({
  ## Check the peculiar function-building code in funxy
  W <- square(1)
  f1a <- function(x, y) sqrt(x^2 + y^2)
  F1a <- funxy(f1a, W)
  f1b <- function(x, y) { sqrt(x^2 + y^2) }
  f2a <- function(x, y) sin(x)
  f2b <- function(x, y) { sin(x) } 
  f3a <- function(x, y) sin(x) + cos(x) 
  f3b <- function(x, y) { sin(x) + cos(x) } 
  f4a <- function(x, y) { z <- x + y ; z }
  f4b <- function(x, y) { x + y }
  F1b <- funxy(f1b, W)
  F2a <- funxy(f2a, W)
  F2b <- funxy(f2b, W)
  F3a <- funxy(f3a, W)
  F3b <- funxy(f3b, W)
  F4a <- funxy(f4a, W)
  F4b <- funxy(f4b, W)
  stopifnot(identical(F1a(cells), F1b(cells)))
  stopifnot(identical(F2a(cells), F2b(cells)))
  stopifnot(identical(F3a(cells), F3b(cells)))
  stopifnot(identical(F4a(cells), F4b(cells)))
  ## check coordinate extraction from objects
  X <- runifrect(9)
  Q <- quadscheme(X)
  a <- F1a(X)
  d <- F1a(Q)
})
}


