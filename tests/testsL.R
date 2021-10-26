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
#'
#'   tests/layered.R
#'
#'   Tests of 'layered' class
#'
#'   $Revision: 1.2 $  $Date: 2020/04/29 08:55:17 $
#'
if(FULLTEST) {
local({
  D <- distmap(cells)
  L <- layered(D, cells,
               plotargs=list(list(ribbon=FALSE), list(pch=16)))
  #'
  plot(L, which=2, plotargs=list(list(pch=3)))
  plot(L, plotargs=list(list(pch=3)))
  #'
  W <- as.owin(L)
  V <- domain(L)
  #' methods
  L2 <- L[square(0.5)]
  Lr <- reflect(L)
  Lf <- flipxy(L)
  Ls <- scalardilate(L, 2)
  La <- shift(L, origin="midpoint")
  Lo <- rotate(L, pi/3, origin="bottomleft")
  Lu <- rescale(L, 0.1, "parsec")
  #' as.layered 
  M <- as.layered(finpines)
  M2 <- as.layered(split(amacrine))
})
}
