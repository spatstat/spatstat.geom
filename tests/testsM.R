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
#'  tests/morpho.R
#' 
#' morphology code blocks
#'
#' $Revision: 1.3 $ $Date: 2020/04/30 02:18:23 $

local({
  if(ALWAYS) { # depends on C code etc
    #' owin
    a <- erosion(letterR, 0.1, polygonal=FALSE)
    b <- dilation(letterR, 0.1, polygonal=FALSE)
    at <- erosion(letterR, 0.1, polygonal=FALSE, strict=TRUE)
    bt <- dilation(letterR, 0.1, polygonal=FALSE, tight=FALSE)
    #' psp
    S <- edges(letterR)
    dm <- dilation(S, 0.1, polygonal=FALSE)
    dt <- dilation(S, 0.1, polygonal=FALSE, tight=FALSE)
    op <- spatstat.options(old.morpho.psp=TRUE)
    dn <- dilation(S, 0.1, polygonal=TRUE)
    spatstat.options(op)
    cS <- closing(S, 0.1, polygonal=FALSE)
    eS <- erosion(S, 0)
    oS <- opening(S, 0)
    #' ppp
    dc <- dilation(cells, 0.06, polygonal=FALSE)
    ec <- erosion(cells, 0)
    oc <- opening(cells, 0)
    #'
    reset.spatstat.options()
  }
})

