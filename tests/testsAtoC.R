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
## badwindowcheck.R
## $Revision: 1.3 $  $Date: 2020/04/28 12:58:26 $
##

local({
  if(ALWAYS) {
    ## Simple example of self-crossing polygon
    x <- read.table("selfcross.txt", header=TRUE)
    ## Auto-repair
    w <- owin(poly=x)

    ## Real data involving various quirks
    b <- read.table("badwindow.txt", header=TRUE)
    b <- split(b, factor(b$i))
    b <- lapply(b, function(z) { as.list(z[,-3]) })
    ## make owin without checking
    W <- owin(poly=b, check=FALSE, fix=FALSE)
    ## Apply stringent checks
    owinpolycheck(W,verbose=FALSE)
    ## Auto-repair
    W2 <- owin(poly=b)
  }
})




##  tests/closeshave.R
## check 'closepairs/crosspairs' code
## validity and memory allocation
## $Revision: 1.28 $ $Date: 2022/02/18 03:06:00 $

## ------- All this code must be run on every hardware -------
local({
  r <- 0.12

  close.all <- closepairs(redwood, r)
  close.ij <- closepairs(redwood, r, what="indices")
  close.ijd <- closepairs(redwood, r, what="ijd")
  close.every <- closepairs(redwood, r, what="all", distinct=FALSE)

  ## test agreement
  stopifnot(identical(close.ij, close.all[c("i","j")]))
  stopifnot(identical(close.ijd, close.all[c("i","j","d")]))

  ## validate basic format of result
  checkformat <- function(object, callstring) {
    if(length(unique(lengths(object))) > 1)
      stop(paste("Result of", callstring,
                 "contains vectors with different lengths"))
    return(invisible(TRUE))
  }
  checkformat(close.all, "closepairs(redwood, r)")  
  checkformat(close.ij, "closepairs(redwood, r, what='indices')")  
  checkformat(close.ijd, "closepairs(redwood, r, what='ijd')")  
  checkformat(close.every,
              "closepairs(redwood, r, what='all', distinct=FALSE)")  

  #' test memory overflow code
  close.cigar <- closepairs(redwood, r, what="ijd", nsize=2)
  close.cigar <- closepairs(redwood, r, what="ijd", nsize=2, periodic=TRUE)

  #' test special cases
  onepoint <- redwood[1]
  checkformat(closepairs(onepoint, r),
              "closepairs(onepoint, r)")
  checkformat(closepairs(onepoint, r, what="indices"),
              "closepairs(onepoint, r, what='indices')")
  checkformat(closepairs(onepoint, r, what="ijd"),
              "closepairs(onepoint, r, what='ijd')")
  checkformat(closepairs(onepoint, r, what="all", distinct=FALSE),
              "closepairs(onepoint, r, what='all', distinct=FALSE)")
  
  #' ..............  crosspairs ..................................
  Y <- split(amacrine)
  on <- Y$on
  off <- Y$off
  
  cross.all <- crosspairs(on, off, r)
  cross.ij <- crosspairs(on, off, r, what="indices")
  cross.ijd <- crosspairs(on, off, r, what="ijd")
  cross.every <- crosspairs(on, off, r, what="all", distinct=FALSE)
  cross.period <- crosspairs(on, off, r, periodic=TRUE)

  ## validate basic format
  checkformat(cross.all, "crosspairs(on, off, r)")
  checkformat(cross.ij, "crosspairs(on, off, r, what='indices')")
  checkformat(cross.ijd, "crosspairs(on, off, r, what='ijd')")
  checkformat(cross.every, "crosspairs(on, off, r, what='all', distinct=FALSE)")
  checkformat(cross.period, "crosspairs(on, off, r, periodic=TRUE)")
  
  ## test agreement
  stopifnot(identical(cross.ij, cross.all[c("i","j")]))
  stopifnot(identical(cross.ijd, cross.all[c("i","j","d")]))

  # closethresh vs closepairs: EXACT agreement
  thresh <- 0.08
  clt <- closethresh(redwood, r, thresh)
  cl <- with(closepairs(redwood, r),
             list(i=i, j=j, th = (d <= thresh)))
  if(!identical(cl, clt))
    stop("closepairs and closethresh disagree")

  reordered <- function(a) {
    o <- with(a, order(i,j))
    as.list(as.data.frame(a)[o,,drop=FALSE])
  }
  samesame <- function(a, b) {
    identical(reordered(a), reordered(b))
  }
  
  ## ...............................................
  #' compare with older, slower code
  op <- spatstat.options(closepairs.newcode=FALSE,
                         closepairs.altcode=FALSE,
                         crosspairs.newcode=FALSE)
  ## ...............................................
  old.close.ij <- closepairs(redwood, r, what="indices")
  old.cross.ij <- crosspairs(on, off, r, what="indices")
  stopifnot(samesame(close.ij, old.close.ij))
  stopifnot(samesame(cross.ij, old.cross.ij))
  # execute only:
  old.close.every <- closepairs(redwood, r, what="all", distinct=FALSE)
  old.close.once <- closepairs(redwood, r, what="all", twice=FALSE)
  #' test memory overflow code
  old.close.cigar <- closepairs(redwood, r, what="ijd", nsize=2)
  old.close.cigar <- closepairs(redwood, r, what="ijd", nsize=2, periodic=TRUE)
  
  ## ...............................................
  spatstat.options(op)
  ## ...............................................

  ## ...............................................
  #' alternative code - execution only
  op <- spatstat.options(closepairs.newcode=FALSE,
                         closepairs.altcode=TRUE)
  alt.close.ij <- closepairs(redwood, r, what="indices")
  alt.close.ijd <- closepairs(redwood, r, what="ijd")
  alt.close.all <- closepairs(redwood, r, what="all")
  #' test memory overflow code
  alt.close.cigar <- closepairs(redwood, r, what="ijd", nsize=2)
  alt.close.cigar <- closepairs(redwood, r, what="ijd", nsize=2, periodic=TRUE)
  spatstat.options(op)
  ## ...............................................
  
  # Rasmus' example
  R <- 0.04
  U <- as.ppp(gridcenters(owin(), 50, 50), W=owin())
  cp <- crosspairs(U, U, R)
  G <- matrix(0, npoints(U), npoints(U))
  G[cbind(cp$i, cp$j)] <- 1
  if(!isSymmetric(G))
    stop("crosspairs is not symmetric in Rasmus example")

  #' periodic distance
  pclose <- function(X, R, method=c("raw", "C")) {
    method <- match.arg(method)
    switch(method,
           raw = {
             D <- pairdist(X, periodic=TRUE)
             diag(D) <- Inf
             result <- which(D <= R, arr.ind=TRUE)
           },
           C = {
             result <- closepairs(X, R, periodic=TRUE, what="indices")
           })
    result <- as.data.frame(result)
    colnames(result) <- c("i","j")
    return(result)
  }
  #' pick a threshold value which avoids GCC bug 323
  RR <- 0.193
  A <- pclose(redwood, RR, "raw")
  B <- pclose(redwood, RR, "C")
  if(!samesame(A,B))
    stop("closepairs.ppp(periodic=TRUE) gives wrong answer")

  #' other functions that don't have a help file
  niets <- crosspairquad(quadscheme(cells), 0.1)

  #' other code blocks
  u <- closepairs(cells, 0.09, periodic=TRUE, what="all")
  v <- closepairs(cells, 0.07, twice=FALSE, neat=TRUE)
  #' tight cluster - guess count does not work
  Xc <- runifrect(100, square(0.01))
  Window(Xc) <- square(1)
  z <- closepairs(Xc, 0.02, what="indices", distinct=FALSE)
  z <- closepairs(Xc, 0.02, what="ijd",     distinct=FALSE)
  z <- closepairs(Xc, 0.02, what="all",     distinct=FALSE)
  #' same task, older code
  aop <- spatstat.options(closepairs.newcode=FALSE)
  z <- closepairs(Xc, 0.02, what="indices", distinct=FALSE)
  z <- closepairs(Xc, 0.02, what="ijd",     distinct=FALSE)
  z <- closepairs(Xc, 0.02, what="all",     distinct=FALSE)
  spatstat.options(aop)

  #' experimental
  r <- 0.08
  a <- closepairs(redwood, r)
  b <- tweak.closepairs(a, r, 26, 0.1, 0.1)
})

local({
  #' Three-dimensional
  ##          X <- runifpoint3(100)
  X <- pp3(runif(100), runif(100), runif(100), box3(c(0,1)))
  cl <- closepairs(X, 0.2, what="indices")
  cl <- closepairs(X, 0.2, what="ijd")
  cl <- closepairs(X, 0.2, distinct=FALSE)
  cl <- closepairs(X, 0.2, distinct=FALSE, what="indices")
  cl <- closepairs(X, 0.2, distinct=FALSE, what="ijd")
  cl <- closepairs(X, 0.2, twice=FALSE, neat=TRUE)
  #' Test memory overflow code
  cl <- closepairs(X, 0.2, what="ijd", nsize=2)
  #' trap obsolete usage
  cl <- closepairs(X, 0.2, ordered=FALSE)
  #' crosspairs
  ##      Y <- runifpoint3(100)
  Y <- pp3(runif(100), runif(100), runif(100), box3(c(0,1)))
  cr <- crosspairs(X, Y, 0.2, what="indices")
  cr <- crosspairs(X, Y, 0.2, what="ijd")
  #' Test memory overflow code
  cr <- crosspairs(X, Y, 0.2, what="ijd", nsize=2)
  #' experimental
  rr <- 0.2
  cl <- closepairs(X, rr)
  ii <- cl$i[[1]]
  xl <- tweak.closepairs(cl, rr, ii, 0.05, -0.05, 0.05)
})

reset.spatstat.options()
#'
#'   tests/cluck.R
#'
#'   Tests of "click*" functions
#'   using queueing feature of spatstatLocator
#'
#'   $Revision: 1.7 $ $Date: 2020/11/02 06:53:30 $

local({
})
## tests/colour.R
##
##  Colour value manipulation and colour maps
##
## $Revision: 1.8 $ $Date: 2020/04/28 08:27:38 $
##

local({
  if(FULLTEST) {
    f <- function(n) grey(seq(0,1,length=n))
    z <- to.grey(f)

    h <- colourmap(rainbow(9), range=c(0.01, 0.1))
    plot(h, labelmap=100)
  }

  if(ALWAYS) {
    a <- colourmap(rainbow(12), range=as.Date(c("2018-01-01", "2018-12-31")))
    print(a)
    print(summary(a))
    a(as.Date("2018-06-15"))

    g <- colourmap(rainbow(4),
                   breaks=as.Date(c("2018-01-01", "2018-04-01",
                                    "2018-07-01", "2018-10-01", "2018-12-31")))
    print(g)
    print(summary(g))
    g(as.Date("2018-06-15"))
  }

  if(FULLTEST) {
    b <- colourmap(rainbow(12), inputs=month.name)
    print(b)
    print(summary(b))
    to.grey(b)
    to.grey(b, transparent=TRUE)
    plot(b, vertical=FALSE)
    plot(b, vertical=TRUE)
    plot(b, vertical=FALSE, gap=0)
    plot(b, vertical=TRUE, gap=0)
    plot(b, vertical=FALSE, xlim=c(0, 2))
    plot(b, vertical=TRUE, xlim=c(0,2))
    plot(b, vertical=FALSE, ylim=c(0, 2))
    plot(b, vertical=TRUE, ylim=c(0,2))

    argh <- list(a="iets", e="niets", col=b, f=42)
    arr <- col.args.to.grey(argh)
    rrgh <- col.args.to.grey(argh, transparent=TRUE)
  }

  if(ALWAYS) {
    #' constant colour map
    colourmap("grey", range=c(0.01, 0.1))
    colourmap("grey", range=as.Date(c("2018-01-01", "2018-12-31")))
    colourmap("grey",
              breaks=as.Date(c("2018-01-01", "2018-04-01",
                               "2018-07-01", "2018-10-01", "2018-12-31")))
    colourmap("grey", inputs=month.name)
  }

  if(FULLTEST) {
    #' empty colour map
    niets <- lut()
    print(niets)
    summary(niets)
    niets <- colourmap()
    print(niets)
    summary(niets)
    plot(niets)
  }

  if(ALWAYS) {
    #' interpolation - of transparent colours
    co <- colourmap(inputs=c(0, 0.5, 1),
                    rgb(red=c(1,0,0), green=c(0,1,0), blue=c(0,0,1),
                        alpha=c(0.3, 0.6, 0.9)))
    plot(interp.colourmap(co))
  }
})

# tests/correctC.R
# check for agreement between C and interpreted code
# for interpoint distances etc.
# $Revision: 1.8 $ $Date: 2020/12/03 03:06:04 $

if(ALWAYS) { # depends on hardware
local({
  eps <- .Machine$double.eps * 4

  checkagree <- function(A, B, blurb) {
    maxerr <- max(abs(A-B))
    cat("Discrepancy", maxerr, "for", blurb, fill=TRUE)
    if(maxerr > eps) 
      stop(paste("Algorithms for", blurb, "disagree"))
    return(TRUE)
  }

  ## pairdist.ppp
  set.seed(190901)
  ## X <- rpoispp(42)
  X <- runifrect(max(2, rpois(1, 42)))
  dC <- pairdist(X, method="C")
  dR <- pairdist(X, method="interpreted")
  checkagree(dC, dR, "pairdist()")

  dCp <- pairdist(X, periodic=TRUE, method="C")
  dRp <- pairdist(X, periodic=TRUE, method="interpreted")
  checkagree(dCp, dRp, "pairdist(periodic=TRUE)")

  dCp2 <- pairdist(X, periodic=TRUE, squared=TRUE, method="C")
  dRp2 <- pairdist(X, periodic=TRUE, squared=TRUE, method="interpreted")
  checkagree(dCp2, dRp2, "pairdist(periodic=TRUE, squared=TRUE)")

  ## crossdist.ppp
  ## Y <- rpoispp(42)
  Y <- runifrect(max(2, rpois(1, 42)))
  dC <- crossdist(X, Y, method="C")
  dR <- crossdist(X, Y, method="interpreted")
  checkagree(dC, dR, "crossdist()")

  dC <- crossdist(X, Y, periodic=TRUE, method="C")
  dR <- crossdist(X, Y, periodic=TRUE, method="interpreted")
  checkagree(dC, dR, "crossdist(periodic=TRUE)")

  dC2 <- crossdist(X, Y, periodic=TRUE, squared=TRUE, method="C")
  dR2 <- crossdist(X, Y, periodic=TRUE, squared=TRUE, method="interpreted")
  checkagree(dC2, dR2, "crossdist(periodic=TRUE, squared=TRUE)")

  # nndist.ppp
  nnC <- nndist(X, method="C")
  nnI <- nndist(X, method="interpreted")
  checkagree(nnC, nnI, "nndist()")

  nn3C <- nndist(X, k=3, method="C")
  nn3I <- nndist(X, k=3, method="interpreted")
  checkagree(nn3C, nn3I, "nndist(k=3)")

  # nnwhich.ppp
  nwC <- nnwhich(X, method="C")
  nwI <- nnwhich(X, method="interpreted")
  checkagree(nwC, nwI, "nnwhich()")

  nw3C <- nnwhich(X, k=3, method="C")
  nw3I <- nnwhich(X, k=3, method="interpreted")
  checkagree(nw3C, nw3I, "nnwhich(k=3)")

  # whist
  set.seed(98123)
  x <- runif(1000)
  w <- sample(1:5, 1000, replace=TRUE)
  b <- seq(0,1,length=101)
  op <- spatstat.options(Cwhist=TRUE)
  aT <- whist(x,b,w)
  spatstat.options(Cwhist=FALSE)
  aF <- whist(x,b,w)
  if(!all(aT == aF))
    stop("Algorithms for whist disagree")
  spatstat.options(op)
})

reset.spatstat.options()
}
