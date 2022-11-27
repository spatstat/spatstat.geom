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
#'   tests/perspim.R
#'
#'   Check persp.im handling of NA, etc
#' 
#'   $Revision: 1.3 $  $Date: 2020/12/04 04:05:54 $

if(FULLTEST) {
local({
  set.seed(42)
  Z <- distmap(letterR, invert=TRUE)[letterR, drop=FALSE]
  X <- runifrect(100, Frame(Z))
  M <- persp(Z, colin=Z, visible=TRUE, phi=50)
  perspPoints(X, Z=Z, M=M)
  P <- psp(c(2.360, 3.079, 2.211),
           c(0.934, 1.881, 2.184),
           c(2.337, 3.654, 3.274),
           c(1.829, 0.883, 2.093), window=letterR)
  perspSegments(P, Z=Z, M=M)
  
  persp(Z, colmap=rainbow)
  persp(Z, colmap=beachcolours, sealevel=mean(Z))
  persp(Z, colin=as.im(Z, dimyx=dim(Z)/4))
})
}
##
## tests/pixelgripes.R
##     Problems related to pixellation of windows
##
## $Revision: 1.8 $ $Date: 2022/10/23 06:21:10 $

if(FULLTEST) {
  local({
    

  ## pixellate.ppp includes mapping from (x,y) to (row, col)
  Z <- pixellate(cells, savemap=TRUE)
  ind <- attr(Z, "map")
  m <- (as.matrix(Z))[ind]
  if(!all(m == 1)) stop("Coordinate mismatch in pixellate.ppp")
})
}
## 
## tests/polygons.R
##
##  $Revision: 1.5 $ $Date: 2020/04/30 05:23:52 $
##
if(ALWAYS) { # involves C code
local({
  co <- as.ppp(corners(letterR), letterR, check=FALSE)
  co[letterR]

  b <- letterR$bdry
  a <- sapply(b, xypolyselfint, yesorno=TRUE)
  a <- lapply(b, xypolyselfint, proper=TRUE)
  
  ## Simple example of self-crossing polygon
  x <- read.table("selfcross.txt", header=TRUE)
  y <- xypolyselfint(x)
})
}
