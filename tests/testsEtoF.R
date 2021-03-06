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
#   $Revision: 1.8 $   $Date: 2020/12/03 03:28:44 $
#
#  Tests of 'funxy' infrastructure etc

if(FULLTEST) {
local({
  ## Check the peculiar function-building code in funxy
  W <- square(1)
  f1a <- function(x, y) sqrt(x^2 + y^2)
  F1a <- funxy(f1a, W)
})
}


