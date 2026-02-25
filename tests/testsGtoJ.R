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
#  $Revision: 1.16 $  $Date: 2025/08/02 08:41:50 $
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

  #' [[ and [[<-
  junk <- pyramidal
  a <- junk[["group"]]
  junk[["group"]] <- sample(a)
  a <- junk[[2]]
  a <- junk[[15,2]]
  junk[[15,2]] <- "schizoaffective"
  junk[[15,2]] <- "z" # Warning given.
  a <- junk[[2]] # The warned-about NA appears as entry 15.
  junk[[10,1]] <- cells
  a <- junk[[10,1]]
  a <- junk[[10,"Neurons"]] 
})
}

if(FULLTEST) {
  local({
    ## NA entries in hyperframes and solists
    h <- hyperframe(A=1:4,
                    B=rep(list(cells), 4),
                    C=letters[1:4],
                    D=rep(list(redwood), 4))
    ## check that column entries all have the same class
    ck <- function(context, x, cls) {
      xname <- deparse(substitute(x))
      if(!all(sapply(x, inherits, what=cls)))
        stop(paste(paste0(context, ","),
                   "some entries of", xname,
                   "do not belong to class", sQuote(cls)))
      invisible(NULL)
    }
    ## cases of "$<-.hyperframe" 
    h$A[2] <- NA
    h$A[2:3] <- NA
    h$B[[2]] <- NA
    ck("After h$B[[2]] <- NA", h$B, "ppp")
    h$B[[2]] <- cells
    ck("After h$B[[2]] <- NA and h$B[[2]] <- cells", h$B, "ppp")
    h$B[2:3] <- NA
    ck("After h$B[2:3] <- NA", h$B, "ppp")
    ## cases of "[<-.hyperframe"
    h[2,"A"] <- NA
    h[2:3, "A"] <- NA
    h[2,"B"] <- NA
    ck("After h[2, 'B'] <- NA", h$B, "ppp")
    h[2,"B"] <- cells
    ck("After h[2, 'B'] <- NA and <- cells", h$B, "ppp")
    h[2:3,"B"] <- NA
    ck("After h[2:3, 'B'] <- NA", h$B, "ppp")
    h[3,] <- NA
    ck("After h[3, ] <- NA", h$B, "ppp")
    stopifnot(all(sapply(h$B, is.ppp)))
    ## cases of solist() and "[<-.solist" and "[[<-.solist"
    x <- solist(cells, NA, cells)
    x[[1]] <- NA
    stopifnot(all(sapply(x, is.ppp)))
    y <- solist(NA, NA, cells)
    stopifnot(all(sapply(y, is.ppp)))
    DC <- as.im(function(x,y){x}, owin())
    z <- solist(DC, DC, NA)
    stopifnot(all(sapply(z, is.im)))
    u <- solist(NA, DC, DC)
    stopifnot(all(sapply(u, is.im)))
    u[[2]] <- NA
    ck("After u[[2]] <- NA", u, "im")
    ## 
    h <- hyperframe(A=1:4,
                    B=rep(list(cells), 4),
                    C=letters[1:4],
                    D=rep(list(redwood), 4))
    shouldfail <- function(words, expr) {
      a <- try(eval(expr), silent=TRUE)
      if(!inherits(a, "try-error"))
        stop(paste(words, "did not generate an error"), call.=FALSE)
      invisible(NULL)
    }
    shouldfail("Assigning an owin into a ppp column",
               {h$B[[2]] <- owin()})
    shouldfail("Assigning two owins into two entries in a ppp column",
               {h$B[2:3] <- rep(list(owin()), 2)})
    
  })
}

if(FULLTEST) {
  #' test internal format of hyperframes in spatstat.data
  #' is consistent with output of 'hyperframe()'
  rebuild <- function(h) {
    rn <- row.names(h)
    do.call(hyperframe, append(as.list(h), list(row.names=rn)))
  }
  stopifnot(identical(cetaceans, rebuild(cetaceans)))
  stopifnot(identical(demohyper, rebuild(demohyper)))
  stopifnot(identical(flu, rebuild(flu)))
  stopifnot(identical(Kovesi, rebuild(Kovesi)))
  stopifnot(identical(osteo, rebuild(osteo)))
  stopifnot(identical(pyramidal, rebuild(pyramidal)))
  stopifnot(identical(simba, rebuild(simba)))
}
#
#  tests/imageops.R
#
#   $Revision: 1.47 $   $Date: 2026/02/24 01:26:07 $
#


if(ALWAYS) {
  local({
    #' Test of case 'CONNECT=24' in src/distmapbin.[ch]
    #' Distance transform with 24-connected neighbours
    A <- distmap(heather$coarse, connect=24)
  })
}

if(FULLTEST) {
local({
  #' cases of 'im' data
  tab <- table(sample(factor(letters[1:10]), 30, replace=TRUE))
  b <- im(tab, xrange=c(0,1), yrange=c(0,10))
  b <- update(b)

  mat <- matrix(sample(0:4, 12, replace=TRUE), 3, 4)
  a <- im(mat)
  levels(a$v) <- 0:4
  a <- update(a)
  
  levels(mat) <- 0:4
  b <- im(mat)
  b <- update(b)

  D <- as.im(mat, letterR)
  df <- as.data.frame(D)
  DD <- as.im(df, step=c(D$xstep, D$ystep))
  
  #' various manipulations
  AA <- A <- as.im(owin())
  BB <- B <- as.im(owin(c(1.1, 1.9), c(0,1)))
  Z <- imcov(A, B)
  stopifnot(abs(max(Z) - 0.8) < 0.1)

  Frame(AA) <- Frame(B)
  Frame(BB) <- Frame(A)
  
  ## handling images with 1 row or column
  
  ycov <- function(x, y) y
  E <- as.im(ycov, owin(), dimyx = c(2,1))
  G <- cut(E, 2)
  H <- as.tess(G)

  E12 <- as.im(ycov, owin(), dimyx = c(1,2))
  G12 <- cut(E12, 2)
  H12 <- as.tess(G12)

  AAA <- as.array(AA)
  EEE <- as.array(E)
  AAD <- as.double(AA)
  EED <- as.double(E)
  aaa <- xtfrm(AAA)
  eee <- xtfrm(E)
  
  ## trivial case in connected.owin
  ## 1 x 1 image
  lego <- im(matrix(1), xrange = c(0, 1), yrange = c(0, 1))
  piecesImage <- connected(lego)
  piecesTess  <- connected(lego, polygonal=TRUE)
  ## cases of 'method' in 'connected.owin'
  d <- distmap(cells, dimyx=32)
  D6 <- (d <= 0.06)
  Z <- connected(D6, method="interpreted")
  Z <- connected(D6,                       connect=4)
  Z <- connected(D6, method="interpreted", connect=4)

  a <- where.max(d, first=FALSE)
  a <- where.min(d, first=FALSE)

  dx <- raster.x(d)
  dy <- raster.y(d)
  dxy <- raster.xy(d)
  xyZ <- raster.xy(Z, drop=TRUE)

  horosho <- conform.imagelist(cells, list(d, Z))

  #' split.im
  W <- square(1)
  X <- as.im(function(x,y){x}, W)
  Y <- dirichlet(runifrect(7, W))
  Z <- split(X, as.im(Y))
  
  ## ...........  cases of "[.im" ........................
  ## index window has zero overlap area with image window
  Out <- owin(c(-0.5, 0), c(0,1))
  oo <- X[Out]
  oo <- X[Out, drop=FALSE]
  if(!is.im(oo)) stop("Wrong format in [.im with disjoint index window")
  oon <- X[Out, drop=TRUE, rescue=FALSE]
  if(is.im(oon)) stop("Expected a vector of values, not an image, from [.im")
  if(!all(is.na(oon))) stop("Expected a vector of NA values in [.im")
  ## 
  Empty <- cells[FALSE]
  ff <- d[Empty]
  gg <- d[2,]
  gg <- d[,2]
  gg <- d[j=2]
  gg <- d[2:4, 3:5]
  hh <- d[2:4, 3:5, rescue=TRUE]
  if(!is.im(hh)) stop("rectangle was not rescued in [.im")
  ## factor and NA values
  f <- cut(d, breaks=4)
  f <- f[f != levels(f)[1], drop=FALSE]
  fff <- f[, , drop=FALSE]
  fff <- f[cells]
  fff <- f[cells, drop=FALSE]
  fff <- f[Empty]
  
  ## ...........  cases of "[<-.im"  .......................
  d[,] <- d[] + 1
  d[Empty] <- 42
  ## smudge() and rasterfilter()
  dd <- smudge(d)

  ## rgb/hsv options
  X <- setcov(owin())
  M <- Window(X)
  Y <- as.im(function(x,y) x, W=M)
  Z <- as.im(function(x,y) y, W=M)
  # convert after rescaling
  RGBscal <- rgbim(X, Y, Z, autoscale=TRUE, maxColorValue=1)
  HSVscal <- hsvim(X, Y, Z, autoscale=TRUE)

  #' cases of [.im
  Ma <- as.mask(M, dimyx=37)
  ZM <- Z[raster=Ma, drop=FALSE]
  ZM[solutionset(Y+Z > 0.4)] <- NA
  ZF <- cut(ZM, breaks=5)
  ZL <- (ZM > 0)
  P <- list(x=c(0.511, 0.774, 0.633, 0.248, 0.798),
            y=c(0.791, 0.608, 0.337, 0.613, 0.819))
  zmp <- ZM[P, drop=TRUE]
  zfp <- ZF[P, drop=TRUE]
  zlp <- ZL[P, drop=TRUE]
  P <- as.ppp(P, owin())
  zmp <- ZM[P, drop=TRUE]
  zfp <- ZF[P, drop=TRUE]
  zlp <- ZL[P, drop=TRUE]

  #' miscellaneous
  ZZ <- zapsmall.im(Z, digits=6)
  ZZ <- zapsmall.im(Z)

  ZS <- shift(Z, origin="centroid")
  ZS <- shift(Z, origin="bottomleft")

  ZA <- affine(Z, mat=diag(c(-1,-2)))

  U <- scaletointerval(Z)
  C <- as.im(1, W=U)
  U <- scaletointerval(C)
  
  #' hist.im
  h <- hist(Z)
  h <- hist(Z, probability=TRUE)
  h <- hist(Z, plot=FALSE)
  Zcut <- cut(Z, breaks=5)
  h <- hist(Zcut) # barplot
  hp <- hist(Zcut, probability=TRUE) # barplot
  plot(h) # plot.barplotdata

  #' plot.im code blocks
  plot(Z, ribside="left")
  plot(Z, ribside="top")
  plot(Z, riblab="value")
  plot(Z, clipwin=square(0.5))
  plot(Z - mean(Z), log=TRUE)
  plot(Z, valuesAreColours=TRUE) # rejected with a warning
  IX <- as.im(function(x,y) { as.integer(round(3*x)) }, square(1))
  co <- colourmap(rainbow(4), inputs=0:3)
  plot(IX, col=co)
  CX <- eval.im(col2hex(IX+1L))
  plot(CX, valuesAreColours=TRUE)
  plot(CX, valuesAreColours=FALSE)

  #' plot.im contour code logarithmic case
  V0 <- setcov(owin())
  V2 <- exp(2*V0+1)
  plot(V2, log=TRUE, addcontour=TRUE, contourargs=list(col="white"))
  plot(V2, log=TRUE, addcontour=TRUE, contourargs=list(col="white", nlevels=2))
  plot(V2, log=TRUE, addcontour=TRUE, contourargs=list(col="white", nlevels=20))
  V4 <- exp(4*V0+1)
  plot(V4, log=TRUE, addcontour=TRUE, contourargs=list(col="white"))
  plot(V4, log=TRUE, addcontour=TRUE, contourargs=list(col="white", nlevels=2))
  plot(V4, log=TRUE, addcontour=TRUE, contourargs=list(col="white", nlevels=20))

  #' pairs.im 
  pairs(solist(Z))
  pairs(solist(A=Z))

  #' harmonise.im (cases)
  H <- harmonise(A=Z)
  H <- harmonise(A=Z, B=Z)
  Zcoarse <- as.im(Z, dimyx=16)
  H <- harmonise(A=Z, B=Zcoarse, D=Z)
  
  #' handling and plotting of character and factor images
  Afactor    <- as.im(col2hex("green"), letterR, na.replace=col2hex("blue"))
  Acharacter <- as.im(col2hex("green"), letterR, na.replace=col2hex("blue"),
                      stringsAsFactors=FALSE)
  plot(Afactor)
  plot(Acharacter, valuesAreColours=TRUE)
  print(summary(Afactor))
  print(summary(Acharacter))

  #' substitute for runifpoint
  rup <- function(n, W) { runifrect(n, Frame(W))[W] }
  #' safelookup (including extrapolation case)
  Z <- as.im(function(x,y) { x - y }, letterR)
  Zcut <- cut(Z, breaks=4)
  B <- grow.rectangle(Frame(letterR), 1)
  X <- superimpose(rup(10, letterR),
                   rup(20, setminus.owin(B, letterR)),
                   vertices(Frame(B)),
                   W=B)
  a <- safelookup(Z, X)
  aa <- safelookup(Z, X, factor=100)
  b <- safelookup(Zcut, X)
  bb <- safelookup(Zcut, X, factor=100)
  cc <- lookup.im(Z, X)
  
  #' im.apply
  #' image data with some NA's
  BE <- bei.extra
  BE[[1]][square(c(200, 300))] <- NA
  #' 'parallel' code: process data as one huge matrix
  Z <- im.apply(BE, sd)  # recognised operation => special code
  Y <- im.apply(BE, max) # unrecognised op => generic code
  #' 'big data' code: break data into several matrices
  spatstat.options(maxmatrix=1200)
  Z2 <- im.apply(BE, sd)
  Y2 <- im.apply(BE, max)
  reset.spatstat.options()
  #' check agreement
  if(max(abs(Y-Y2)) > sqrt(.Machine$double.eps))
    stop("Inconsistent results from im.apply max (parallel vs batch-wise)")
  if(max(abs(Z-Z2)) > sqrt(.Machine$double.eps))
    stop("Inconsistent results from im.apply sd (parallel vs batch-wise)")
  #' check unitname information was retained
  metres <- unitname(BE[[1L]])  
  chuk <- function(A, txt) {
    if(!identical(unitname(A), metres))
      stop(paste(txt, "does not retain unitname"))
  }
  chuk(Z, "im.apply(bei.extra, sd)")
  chuk(Y, "im.apply(bei.extra, max)")
  chuk(Z2, "im.apply(bei.extra, sd) for big data")
  chuk(Y2, "im.apply(bei.extra, max) for big data")
  
  #' Math.imlist, Ops.imlist, Complex.imlist
  U <- Z+2i
  B <- U * (2+1i)
  print(summary(B))
  V <- solist(A=U, B=B)
  negV <- -V
  E <- Re(V)
  negE <- -E

})
}

if(ALWAYS) {
  local({
    #' check nearest.valid.pixel
    W <- Window(demopat)
    set.seed(911911)
    X <- runifrect(1000, Frame(W))[W]
    Z <- quantess(W, function(x,y) { x }, 9)$image
    nearest.valid.pixel(numeric(0), numeric(0), Z)
    x <- X$x
    y <- X$y
    a <- nearest.valid.pixel(x, y, Z, method="interpreted")
    b <- nearest.valid.pixel(x, y, Z, method="C")
    if(!isTRUE(all.equal(a,b)))
      stop("Unequal results in nearest.valid.pixel")
      if(!identical(a,b)) 
        stop("Equal, but not identical, results in nearest.valid.pixel")
  })
}

if(ALWAYS) {
  local({
    #' rastersample
    Z <- distmap(letterR, invert=TRUE, eps=0.02) ## non-square raster
    Zs <- rastersample(Z, phase=2, scale=3)
    Zf <- rastersample(Z, phase=1.5, scale=2.7)
    Zss <- rastersample(Z, Zs)
    Zff <- rastersample(Z, Zf)
  })
}
    


