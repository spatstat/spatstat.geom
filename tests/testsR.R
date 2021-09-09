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
#'  tests/randbasic.R
#'   Tests of basic random generation code 
#'  $Revision: 1.1 $ $Date: 2021/09/09 09:59:23 $


local({
  if(FULLTEST) {
    #' cases not covered in examples
    A <- runifrect(6, nsim=2)
    A <- rsyst(nx=4, nsim=2)
    A <- rjitter(cells, nsim=2, retry=FALSE)
    A <- rjitter(cells, nndist(cells)/2)
    A <- rjitter(cells[FALSE])
  }
})


