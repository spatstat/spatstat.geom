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
#    $Revision: 1.5 $  $Date: 2020/04/28 12:58:26 $
#

if(ALWAYS) {
local({
  ## make a factor image
  m <- factor(rep(letters[1:4], 4))
  Z <- im(m, xcol=1:4, yrow=1:4)
  ## make a point pattern
  set.seed(42)
  X <- runifpoint(20, win=as.owin(Z))
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
#' tests/formuli.R
#'
#'  Test machinery for manipulating formulae
#' 
#' $Revision: 1.7 $  $Date: 2020/04/28 12:58:26 $

local({

  ff <- function(A, deletevar, B) {
    D <- reduceformula(A, deletevar)
    if(!spatstat.utils::identical.formulae(D, B)) {
      AD <- as.expression(substitute(reduceformula(A,d),
                                     list(A=A, d=deletevar)))
      stop(paste(AD, "\n\tyields ", spatstat.utils::pasteFormula(D),
                 " instead of ", spatstat.utils::pasteFormula(B)),
           call.=FALSE)
    }
    invisible(NULL)
  }

  ff(~ x + z, "x", ~z)

  ff(y ~ x + z, "x", y~z)

  ff(~ I(x^2) + z, "x",  ~z)

  ff(y ~ poly(x,2) + poly(z,3), "x", y ~poly(z,3))

  ff(y ~ x + z, "g", y ~ x + z)

  reduceformula(y ~ x+z, "g", verbose=TRUE)
  reduceformula(y ~ sin(x-z), "z", verbose=TRUE)
  
  illegal.iformula(~str*g, itags="str", dfvarnames=c("marks", "g", "x", "y"))
})



#
#  tests/func.R
#
#   $Revision: 1.7 $   $Date: 2020/11/02 07:00:08 $
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


