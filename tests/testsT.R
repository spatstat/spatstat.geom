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
#'   tests/tessera.R
#'   Tessellation code, not elsewhere tested
#'   $Revision: 1.9 $ $Date: 2020/12/04 08:04:38 $
#'
if(FULLTEST) {
local({
  W <- owin()
  Wsub <- square(0.5)
  X <- runifrect(7, W)
  A <- dirichlet(X)
  marks(A) <- 1:nobjects(A)
  Z <- distmap(letterR, invert=TRUE)[letterR, drop=FALSE]
  H <- tess(xgrid=0:2, ygrid=0:3)
  #' discretisation of tiles
  V <- as.im(A)
  B <- tess(window=as.mask(W), tiles=tiles(A))
  #' logical images
  D <- tess(image=(Z > 0.2))
  U <- (Z > -0.2) # TRUE or NA
  E <- tess(image=U, keepempty=TRUE)
  G <- tess(image=U, keepempty=FALSE)
  #' methods
  flay <- function(op, ..., Rect=H, Poly=A, Img=E) {
    a <- do.call(op, list(Rect, ...))
    b <- do.call(op, list(Poly, ...))
    e <- do.call(op, list(Img, ...))
  }
  flay(reflect)
  flay(flipxy)
  flay(shift, vec=c(1,2))
  flay(scalardilate, f=2) 
  flay(rotate, angle=pi/3, centre=c(0, 0))
  flay(rotate, angle=pi/2)
  flay(affine, mat=matrix(c(1,2,0,1), 2, 2), vec=c(1,2))
  flay(affine, mat=diag(c(1,2)))
  flay(as.data.frame)
  ##
  unitname(A) <- "km"
  unitname(B) <- c("metre", "metres")
  unitname(B)
  print(B)
  Bsub <- B[c(3,5,7)]
  print(Bsub)
  tilenames(H) <- letters[seq_along(tilenames(H))]
  G <- tess(xgrid=(0:3)/3, ygrid=(0:3)/3)
  tilenames(G) <- letters[1:9]
  h <- tilenames(G)
  GG <- as.tess(tiles(G))
  #'
  Pe <- intersect.tess(A, Wsub, keepmarks=TRUE)
  Pm <- intersect.tess(A, as.mask(Wsub), keepmarks=TRUE)
  H <- dirichlet(runifrect(4, W))
  AxH <- intersect.tess(A, H, keepmarks=TRUE) # A is marked, H is not
  HxA <- intersect.tess(H, A, keepmarks=TRUE) # A is marked, H is not
  
  b <- bdist.tiles(D)
  b <- bdist.tiles(A[c(3,5,7)])
  #'
  Eim <- as.im(E, W=letterR)
  #'
  #' chop.tess
  #'    horiz/vert lines
  W <- square(1)
  H <- infline(h=(2:4)/5)
  V <- infline(v=(3:4)/5)
  WH <- chop.tess(W, H)
  WV <- chop.tess(W, V)
  #'     polygonal tessellation
  D <- dirichlet(runifrect(4))
  DH <- chop.tess(D, H)
  DV <- chop.tess(D, V)
  #'     image-based tessellation
  f <- function(x,y){factor(round(4* (x^2 + y^2)))}
  A <- tess(image=as.im(f, W=W))
  L <- infline(p=(1:3)/3, theta=pi/4)
  AL <- chop.tess(A, L)
  AH <- chop.tess(A, H)
  AV <- chop.tess(A, V)
  #'
  #' quantess
  #' quantess.owin
  a <- quantess(square(1), "x", 3)
  a <- quantess(square(1), "y", 3)
  a <- quantess(square(1), "rad", 5, origin=c(1/2, 1/3))
  a <- quantess(square(1), "ang", 7, origin=c(1/2, 1/3))
  ZFUN <- function(x,y){y-x}
  a <- quantess(square(1), ZFUN, 3)
  b <- quantess(letterR, "y", 3)
  #' quantess.ppp
  d <- quantess(cells, "y", 4)
  g <- quantess(demopat, "x", 5)
  g <- quantess(demopat, "y", 5)
  g <- quantess(demopat, "rad", 5, origin=c(4442, 4214))
  g <- quantess(demopat, "ang", 5, origin=c(4442, 4214))
  g <- quantess(demopat, ZFUN, 7)
  #' quantess.im
  D <- distmap(demopat)
  h <- quantess(D, "y", 4)
  h <- quantess(D, ZFUN, 5)
  g <- quantess(D, "rad", 5, origin=c(4442, 4214))
  g <- quantess(D, "ang", 5, origin=c(4442, 4214))
  #'
  X <- shift(chorley, vec = c(1e6, 0))
  tes <- quantess(X, "x", 4)
  if(anyDuplicated(tilenames(tes)))
    stop("quantess produced non-unique tilenames")
  ## 
  ##
  XR <- runifrect(40, Frame(letterR))[letterR]
  da <- dirichletAreas(discretise(XR))
})
}
#'    tests/trigraph.R
#'
#'   Tests for C code in trigraf.c
#'   
#'  $Revision: 1.5 $  $Date: 2020/06/12 00:35:44 $
#'
if(ALWAYS) { # depends on C code 
local({
  #' called from deldir.R
  spatstat.deldir.setopt(FALSE, TRUE)
  A <- delaunay(redwood)
  spatstat.deldir.setopt(FALSE, FALSE)
  B <- delaunay(redwood)
  spatstat.deldir.setopt(TRUE, TRUE)
  #' called from edges2triangles.R
  tryangles <- function(iedge, jedge, nt=0) {
    spatstat.options(fast.trigraph=FALSE)
    A <- edges2triangles(iedge, jedge)
    spatstat.options(fast.trigraph=TRUE)
    B <- edges2triangles(iedge, jedge)
    if(!all(dim(A) == dim(B)) || !all(A == B))
      stop(paste("Discrepancy in edges2triangles (with", nt, "triangles)"))
  }
  ## ii <- simplenet$from
  ## jj <- simplenet$to
  ii <- c(1, 3, 4, 2, 4, 5, 5, 6, 7, 8)
  jj <- c(4, 4, 5, 6, 6, 8, 9, 10, 10, 10)
  tryangles(ii,          jj,          0)
  tryangles(c(ii, 1),    c(jj, 5),    1)
  tryangles(c(ii, 1, 8), c(jj, 5, 9), 2)
})
}
reset.spatstat.options()


