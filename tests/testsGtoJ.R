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
#
# tests/hyperframe.R
#
# test "[.hyperframe" etc
#
#  $Revision: 1.9 $  $Date: 2020/12/03 03:32:13 $
#

if(FULLTEST) {
local({
  lambda <- runif(4, min=50, max=100)
  X <- lapply(as.list(lambda), function(x) { runifrect(rpois(1, x)) })
  h <- hyperframe(lambda=lambda, X=X)
  h$lambda2 <- lambda^2
  h[, "lambda3"] <- lambda^3
  h[, "Y"] <- X
  h[, "X"] <- lapply(X, flipxy)
  h[, c("X", "Y")] <- hyperframe(X=X, Y=X)

  names(h) <- LETTERS[1:5]
  print(h)

  summary(h)
  str(h)
  head(h)
  tail(h)

  rn <- rownames(h)
  r.n <- row.names(h)
  if(!identical(rn, r.n))
    stop("rownames and row.names conflict for hyperframes")

  dn <- dimnames(h)
  dimnames(h) <- dn
  dimnames(h)[[2]][2] <- "copacetic"
  dimnames(h)[[1]][2] <- "second"

  #' hyperframe with a hyperatom
  H <- hyperframe(A=runif(3), B=1:3, D=runifrect(10))
  H[,3]
  H[,3,drop=TRUE]
  #' special cases of [<-
  H$B <- H[,1]
  H[2:3,1] <- H[2:3,2]
  H[2:3,1] <- H[2,2]
  H[2,1:2] <- H[3,1:2]

  #' split
  f <- factor(c("a", "a", "b"))
  G <- split(H, f)
  G[["a"]]$B <- 42
  split(H, f) <- G
})
}
#
#  tests/imageops.R
#
#   $Revision: 1.35 $   $Date: 2022/04/14 00:49:39 $
#


