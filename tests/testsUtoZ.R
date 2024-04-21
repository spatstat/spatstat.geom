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
# tests/utils.R
#
# Tests of miscellaneous utilities 
#
#  $Revision: 1.1 $  $Date: 2023/05/07 08:59:32 $

local({
  if(FULLTEST) {
    ## test code blocks in 'progressreport'
    pstate <- list()
    for(i in 1:10) {
      Sys.sleep(1)
      pstate <- progressreport(i, 10,
                               formula = (time ~ i + I(i^2) + I(i^3)),
                               showtime=TRUE,
                               savehistory=TRUE,
                               state=pstate)
    }
  }
})

#
# tests/windows.R
#
# Tests of owin geometry code
#
#  $Revision: 1.18 $  $Date: 2023/07/11 06:21:37 $

local({
  if(ALWAYS) { # C code
    ## Ege Rubak spotted this problem in 1.28-1
    A <- as.owin(ants)
    B <- dilation(A, 140)
    if(!is.subset.owin(A, B))
      stop("is.subset.owin fails in polygonal case")

    ## thanks to Tom Rosenbaum
    A <- shift(square(3), origin="midpoint")
    B <- shift(square(1), origin="midpoint")
    AB <- setminus.owin(A, B)
    D <- shift(square(2), origin="midpoint")
    if(is.subset.owin(D,AB))
      stop("is.subset.owin fails for polygons with holes")

    ## thanks to Brian Ripley / SpatialVx
    M <- as.mask(letterR)
    stopifnot(area(bdry.mask(M)) > 0)
    stopifnot(area(convexhull(M)) > 0)
    R <- as.mask(square(1))
    stopifnot(area(bdry.mask(R)) > 0)
    stopifnot(area(convexhull(R)) > 0)
  }

  if(FULLTEST) {
    RR <- convexify(as.mask(letterR))
    CC <- covering(letterR, 0.05, eps=0.1)
  
    #' as.owin.data.frame
    V <- as.mask(letterR, eps=0.2)
    Vdf <- as.data.frame(V)
    Vnew <- as.owin(Vdf)
    zz <- mask2df(V)
  }

  if(ALWAYS) { # C code
    RM <- owinpoly2mask(letterR, as.mask(Frame(letterR)), check=TRUE)
  }

  if(FULLTEST) {
    #' as.owin
    U <- as.owin(quadscheme(cells))
    U2 <- as.owin(list(xmin=0, xmax=1, ymin=0, ymax=1))
  }

  if(ALWAYS) {
    #' validity of as.mask applied to rectangles with additional raster info
    Z <- as.im(unit.square())
    R <- square(0.5)
    aR <- area(R)
    a <- area(as.mask(R, xy=Z))
    if(abs(a-aR) > aR/20)
      stop("Problem with as.mask(rectangle, xy=image)")
    a <- area(as.mask(R, xy=list(x=Z$xcol, y=Z$yrow)))
    if(abs(a-aR) > aR/20)
      stop("Problem with as.mask(rectangle, xy=list(x,y))")
  }
    
  if(FULLTEST) {
    #' intersections involving masks
    B1 <- square(1)
    B2 <- as.mask(shift(B1, c(0.2, 0.3)))
    o12 <- overlap.owin(B1, B2)
    o21 <- overlap.owin(B2, B1)
    i12 <- intersect.owin(B1, B2, eps=0.01)
    i21 <- intersect.owin(B2, B1, eps=0.01)
    E2 <- emptywindow(square(2))
    e12 <- intersect.owin(B1, E2)
    e21 <- intersect.owin(E2, B1)
  
    #' geometry
    inradius(B1)
    inradius(B2)
    inradius(letterR)
    inpoint(B1)
    inpoint(B2)
    inpoint(letterR)
    is.convex(B1)
    is.convex(B2)
    is.convex(letterR)
    volume(letterR)
    perimeter(as.mask(letterR))
    boundingradius(cells)
  
    boundingbox(letterR)
    boundingbox(letterR, NULL)
    boundingbox(solist(letterR))

  }

  if(ALWAYS) { # C code
    spatstat.options(Cbdrymask=FALSE)
    bb <- bdry.mask(letterR)
    spatstat.options(Cbdrymask=TRUE)
  }

  if(FULLTEST) {
    X <- longleaf[square(50)]
    marks(X) <- marks(X)/8
    D <- discs(X)
    D <- discs(X, delta=5, separate=TRUE)
  }

  if(ALWAYS) { # C code
    AD <- dilated.areas(cells,
                        r=0.01 * matrix(1:10, 10,1),
                        constrained=FALSE, exact=FALSE)
  }

  if(FULLTEST) {
    periodify(B1, 2)
    periodify(union.owin(B1, B2), 2)
    periodify(letterR, 2)
  }

  if(ALWAYS) {
    #' Ancient bug in inside.owin
    W5 <- owin(poly=1e5*cbind(c(-1,1,1,-1),c(-1,-1,1,1)))
    W6 <- owin(poly=1e6*cbind(c(-1,1,1,-1),c(-1,-1,1,1)))
    i5 <- inside.owin(0,0,W5)
    i6 <- inside.owin(0,0,W6)
    if(!i5) stop("Wrong answer from inside.owin")
    if(i5 != i6) stop("Results from inside.owin are scale-dependent")
  }

  if(FULLTEST) {
    #' miscellaneous utilities
    thrash <- function(f) {
      f(letterR)
      f(Frame(letterR))
      f(as.mask(letterR))
    }
    thrash(meanX.owin)
    thrash(meanY.owin)
    thrash(intX.owin)
    thrash(intY.owin)

    interpretAsOrigin("right", letterR)
    interpretAsOrigin("bottom", letterR)
    interpretAsOrigin("bottomright", letterR)
    interpretAsOrigin("topleft", letterR)
    interpretAsOrigin("topright", letterR)
  }

  if(ALWAYS) { # depends on polyclip
    A <- break.holes(letterR)
    B <- break.holes(letterR, splitby="y")
    plot(letterR, col="blue", use.polypath=FALSE)
  }

  if(ALWAYS) {  # C code
    #' mask conversion
    M <- as.mask(letterR)
    D2 <- as.data.frame(M)              # two-column
    D3 <- as.data.frame(M, drop=FALSE)  # three-column
    M2 <- as.owin(D2)
    M3 <- as.owin(D3)
    W2 <- owin(mask=D2)
    W3 <- owin(mask=D3)
  }

  if(FULLTEST) {
    #' void/empty cases
    nix <- nearest.raster.point(numeric(0), numeric(0), M)
    E <- emptywindow(Frame(letterR))
    print(E)
    #' cases of summary.owin
    print(summary(E)) # empty
    print(summary(Window(humberside))) # single polygon
    #' additional cases of owin()
    B <- owin(mask=M$m) # no pixel size or coordinate info
    xy <- as.data.frame(letterR)
    xxyy <- split(xy[,1:2], xy$id)
    spatstat.options(checkpolygons=TRUE)
    H <- owin(poly=xxyy, check=TRUE)
  }

  #' Code for/using intersection and union of windows

  if(FULLTEST) {
    Empty <- emptywindow(Frame(letterR))
    a <- intersect.owin()
    a <- intersect.owin(Empty)
    a <- intersect.owin(Empty, letterR)
    a <- intersect.owin(letterR, Empty)
    b <- intersect.owin()
    b <- intersect.owin(Empty)
    b <- intersect.owin(Empty, letterR)
    b <- intersect.owin(letterR, Empty)
    d <- union.owin(as.mask(square(1)), as.mask(square(2)))
    #' [.owin
    A <- erosion(letterR, 0.2)
    Alogi <- as.im(TRUE, W=A)
    B <- letterR[A]
    B <- letterR[Alogi]
    #' miscellaneous
    D <- convexhull(Alogi)
  }
})

reset.spatstat.options()
##
## tests/xysegment.R
##                      [SEE ALSO tests/segments.R]
##
##    Test weird problems and boundary cases for line segment code
##
##    $Version$ $Date: 2022/10/23 01:21:09 $ 
##

local({
  if(FULLTEST) {
    ## segment of length zero
    B <- psp(1/2, 1/2, 1/2, 1/2, window=square(1))
    BB <- angles.psp(B)
    A <- runifrect(3)
    AB <- project2segment(A,B)

    ## mark inheritance
    X <- psp(runif(10), runif(10), runif(10), runif(10), window=owin())
    marks(X) <- 1:10
    Y <- selfcut.psp(X)
    marks(X) <- data.frame(A=1:10, B=factor(letters[1:10]))
    Z <- selfcut.psp(X)
    #' psp class support
    S <- unmark(X)
    marks(S) <- sample(factor(c("A","B")), nobjects(S), replace=TRUE)
    intensity(S)
    intensity(S, weights=runif(nsegments(S)))
  }
})


