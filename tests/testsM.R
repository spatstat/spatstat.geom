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
#' $Revision: 1.4 $ $Date: 2026/07/07 05:57:51 $

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

  if(FULLTEST) {
    #' example from Alexey Sergushichev
    pixelSize <- 0.1
    #' 3x3 cross inside 40x40 pixel raster window.
    croix <- as.mask(owin(xrange = c(0, 4), yrange = c(0, 4)), eps = pixelSize)
    croix$m[] <- FALSE
    croix$m[33:35, 34] <- TRUE
    croix$m[34, 33:35] <- TRUE
    #' check the cross has 5 pixels
    nblack <- sum(as.im(croix))
    stopifnot(nblack == 5)
    ## try an increasing sequence of small radii 
    rvals <- ((0:10)/5) * pixelSize
    counts <- sapply(rvals, function(r) sum(as.im(dilation(croix, r))))
    print(data.frame(rvals, counts))
    if(any(counts < nblack))
      stop("Small dilation(s) failed to contain input region")
    if(any(diff(counts) < 0))
      stop("Area of dilation was not a monotone function of radius")
  }

})

