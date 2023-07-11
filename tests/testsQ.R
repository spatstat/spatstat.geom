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
#'    tests/quadschemes.R
#'
#'    Quadrature schemes, dummy points etc
#' 
#'   $Revision: 1.8 $ $Date: 2020/12/04 04:56:26 $
#'

if(FULLTEST) {
local({
  ##  class 'quad' 
  qu <- quadscheme(cells)
  qm <- quadscheme(amacrine)
  plot(qu)
  plot(qm)
  is.multitype(qu)
  is.multitype(qm)
  a <- param.quad(qu)
  a <- param.quad(qm)
  a <- equals.quad(qu)
  a <- equals.quad(qm)
  a <- domain(qu)
  unitname(qu) <- c("Furlong", "Furlongs")
  
  ## utilities
  b <- cellmiddles(square(1), 3, 4)
  b <- cellmiddles(letterR, 3, 4, distances=FALSE)
  b <- cellmiddles(letterR, 3, 4, distances=TRUE)
  v <- tilecentroids(square(1), 3, 4)
  v <- tilecentroids(letterR, 3, 4)
  n <- default.n.tiling(cells)
  n <- default.n.tiling(cells, nd=4)
  n <- default.n.tiling(cells, ntile=4)
  n <- default.n.tiling(cells, ntile=4, quasi=TRUE)

  ## quadrature weights - special cases
  ## X <- runifpoint(10, as.mask(letterR))
  X <- runifrect(10, Frame(letterR))[as.mask(letterR)]
  gr <- gridweights(X, ntile=12, npix=7) # causes warnings about zero digital area
  
  ## plot.quad 
  plot(quadscheme(cells, method="dirichlet", nd=7),              tiles=TRUE)
  plot(quadscheme(cells, method="dirichlet", nd=7, exact=FALSE), tiles=TRUE)
  
  ## logistic
  d <- quadscheme.logi(cells, logi.dummy(cells, "binomial"))
  print(summary(d))
  d <- quadscheme.logi(cells, logi.dummy(cells, "poisson"))
  print(summary(d))
  d <- quadscheme.logi(cells, logi.dummy(cells, "grid"))
  print(summary(d))
  d <- quadscheme.logi(cells, logi.dummy(cells, "transgrid"))
  print(summary(d))
  d <- quadscheme.logi(amacrine,
                       logi.dummy(amacrine, "binomial", mark.repeat=TRUE))
  print(summary(d))
  d <- quadscheme.logi(amacrine,
                       logi.dummy(amacrine, "poisson", mark.repeat=FALSE))
  print(summary(d))
})
}
#
# tests/quadcount.R
#
# Tests of quadrat counting code
#
#  $Revision: 1.1 $  $Date: 2023/07/11 06:55:50 $

local({
  if(FULLTEST) {
    ## from M. Gimond
    A <- quadratcount(humberside, 2, 3)
    nA <- as.integer(A)
    if(!all(nA == c(2, 13, 34, 20, 11, 123)))
      stop("Incorrect quadrat count (2,3)")
    ## execute intensity.quadratcount
    lamA <- intensity(A, image=TRUE)
    ## check sum 1/lambda equals area
    vA <- sum(1/lamA[humberside])
    if(abs(1 - vA/area(Window(humberside))) > 0.05)
      stop("Incorrect sum of 1/lambda (2,3)")
    ##
    B <- quadratcount(humberside, 5, 3)
    nB <- as.integer(B) 
    if(!all(nB == c(0, 3, 0, 2, 2, 3, 14, 117, 19, 5, 35, 0, 3)))
      stop("Incorrect quadrat count (5,3)")
    lamB <- intensity(B, image=TRUE)
    vB <- sum(1/lamB[humberside])
    aaB <- tile.areas(as.tess(B))
    tilemap <- attr(B, "tilemap")
    retaintile <- logical(length(aaB))
    retaintile[tilemap] <- (nB > 0)
    aB <- sum(aaB[retaintile])
    if(abs(1 - vB/aB) > 0.05)
      stop("Incorrect sum of 1/lambda (5,3)")
  }
})
reset.spatstat.options()
