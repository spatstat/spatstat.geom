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
#'   tests/ppp.R
#'
#'   $Revision: 1.12 $ $Date: 2020/12/04 04:43:33 $
#'
#'  Untested cases in ppp() or associated code

local({
  ## X <- runifpoint(10, letterR)
  ## Y <- runifpoint(10, complement.owin(letterR))
  Bin  <- owin(c(2.15, 2.45), c(0.85, 3.0))
  Bout <- owin(c(2.75, 2.92), c(0.85, 1.4))
  X <- runifrect(10, Bin)[letterR]
  Y <- runifrect(10, Bout)[complement.owin(letterR)]
  
  if(FULLTEST) {
    #' test handling of points out-of-bounds
    df <- rbind(as.data.frame(X), as.data.frame(Y))
    A <- ppp(df$x, df$y, window=letterR, marks=1:20)
    #' test handling of points with bad coordinates
    df$x[1:3] <- c(Inf, NA, NaN)
    df$y[18:20] <- c(Inf, NA, NaN)
    B <- ppp(df$x, df$y, window=letterR, marks=1:20)
    D <- ppp(df$x, df$y, window=letterR, marks=data.frame(id=1:20, u=runif(20)))
  
    #' test print/summary/plot methods on these bad objects
    print(A)
    print(B)
    print(D)
    print(summary(A))
    print(summary(B))
    print(summary(D))
    plot(A)
    plot(B)
    plot(D)
    plot(attr(A, "rejects"))
    plot(attr(B, "rejects"))
    plot(attr(D, "rejects"))
  
    #' subset operator --- cases not covered elsewhere
    #'   subset index is a logical image
    Z <- distmap(letterR, invert=TRUE)
    V <- (Z > 0.2)
    XV <- X[V]
    #'   multiple columns of marks
    fun3 <- finpines[1:3]
    #'   multiple columns of marks, one of which is a factor
    U <- finpines
    marks(U)[,2] <- factor(c(rep("A", 60), rep("B", npoints(U)-60)))
    UU <- U[1:3, drop=TRUE]

    #' cut.ppp
    CU <- cut(U, "height")
    CU <- cut(U, breaks=3)

    #' cases of [<-.ppp
    set.seed(999)
    X <- cells
    B <- square(0.2)
    X[B] <- runifrect(3, B)
    #' checking 'value'
    Y <- flipxy(X)
    X[B] <- Y[square(0.3)]
    ## deprecated use of second argument
    X[,1:4] <- runifrect(3)  # deprecated
    X[,B] <- runifrect(3, B) # deprecated 
    X[1:3, B] <- runifrect(20)
    A <- superimpose(cells, X, W="convex")
    A <- superimpose(cells, X, W=ripras)
    B <- superimpose(concatxy(cells), concatxy(X), W=NULL)
    ## superimpose.splitppp
    Y <- superimpose(split(amacrine))

    ## catch outdated usage of scanpp
    d <- system.file("rawdata", "amacrine", package="spatstat.data")
    if(nzchar(d)) {
      W <- owin(c(0, 1060/662), c(0, 1))
      Y <- scanpp("amacrine.txt", dir=d, window=W, multitype=TRUE)
      print(Y)
    }
    ## (bad) usage of cobble.xy
    xx <- runif(10)
    yy <- runif(10)
    W1 <- cobble.xy(xx, yy)
    W2 <- cobble.xy(xx, yy, boundingbox)
    Wnope <- cobble.xy(xx, yy, function(x,y) {cbind(x,y)}, fatal=FALSE)
  }
})

#
# tests/ppx.R
#
# Test operations for ppx objects
#
#  $Revision: 1.9 $ $Date: 2020/12/04 04:49:40 $
#

local({
  if(ALWAYS) {
    ## make data
    df <- data.frame(x=c(1,2,2,1)/4, y=c(1,2,3,1)/4, z=c(2,3,4,3)/5)
    X <- ppx(data=df, coord.type=rep("s", 3), domain=box3())
  }
  if(ALWAYS) {
    #' methods involving C code
    unique(X)
    duplicated(X)
    anyDuplicated(X)
    multiplicity(X)
    uniquemap(X)
  }
  if(FULLTEST) {
    #' general tests
    print(X)
    summary(X)
    plot(X)
    domain(X)
    unitname(X) <- c("metre", "metres")
    unitname(X)

    #' subset operator
    X[integer(0)]
    Y <- X %mark% data.frame(a=df$x, b=1:4)
    Y[1:2]
    Y[FALSE]
    marks(Y) <- as.data.frame(marks(Y))
    Y[integer(0)]
    Y[1:2]
    Y[FALSE]
  }

  if(FULLTEST) {
    #' two dimensional
    A <- ppx(data=df[,1:2], coord.type=rep("s", 2), domain=square(1))
    plot(A)
    B <- ppx(data=df[,1:2], coord.type=rep("s", 2), domain=NULL)
    plot(B)
    #' one dimensional
    E <- ppx(data=data.frame(x=runif(10)))
    plot(E)
  
    #' bug
    stopifnot(identical(unmark(chicago[1]),
                        unmark(chicago)[1]))

    #' ppx with zero points
    U <- chicago[integer(0)]
    V <- U %mark% 1
    V <- U %mark% factor("a")

    #' simplify lower-dimensional patterns
    X3 <- ppx(data=df, coord.type=rep("s", 3), domain=box3(), simplify=TRUE)
    stopifnot(is.pp3(X3))
    X2 <- ppx(data=df[,1:2], coord.type=rep("s", 2), domain=square(1), simplify=TRUE)
    stopifnot(is.ppp(X2))

    #' marks<-.ppx
    M <- as.matrix(X)
    marks(X) <- df[,1]
    marks(X) <- df[,integer(0)]
  }

  if(FULLTEST) {
    ## ............  from Ege ..........................
    ## Tests for shift:
    ## Check ppp and ppx shift are the same
    X <- cells
    Y <- ppx(coords(cells), domain = boxx(0:1,0:1))
    Xs <- shift(X, vec = c(1,1))
    Ys <- shift(Y, vec = c(1,1))
    stopifnot(all.equal(coords(Xs), coords(Ys),
                        check.attributes = FALSE))
    stopifnot(all.equal(domain(Xs), as.owin(domain(Ys)),
                        check.attributes = FALSE))
    ## Check a single numeric for vec in shift.ppx
    stopifnot(identical(Ys, shift(Y, vec = 1)))

    ## Tests for scale:
    dat <- data.frame(x=1:3, y=1:3, m=letters[1:3])
    xrange <- yrange <- c(0,4)
    cent <- c(2,2)
    scal <- c(5,5)
    X <- as.ppp(dat, W = owin(xrange, yrange))
    Xscaled <- affine(shift(X, vec = -cent), mat = diag(1/scal))
    ## Check ppx without domain:
    Y <- ppx(dat, coord.type = c("spatial", "spatial", "mark"))
    Yscaled <- scale(Y, center = cent, scale = scal)
    stopifnot(all.equal(coords(Xscaled),
                        coords(Yscaled),
                        check.attributes = FALSE))
    ## Check ppx with domain:
    Y$domain <- boxx(xrange, yrange)
    Yscaled <- scale(Y, center = cent, scale = scal)
    stopifnot(all.equal(as.boxx(Window(Xscaled)),
                        domain(Yscaled),
                        check.attributes = FALSE))

    ## Tests for intersect.boxx:
    ## Should be unit 2D box:
    A <- intersect.boxx(boxx(c(-1,1),c(0,2)), boxx(c(0,3),c(0,1)))
    stopifnot(identical(A, boxx(c(0,1),c(0,1))))
    ## Should be empty (NULL)
    B <- intersect.boxx(boxx(c(-1,1),c(0,2)),
                        boxx(c(0,3),c(0,1)),
                        boxx(c(1,2), c(-1,1)))
    stopifnot(is.null(B))
    ## Should be unit 3D box:
    C <- intersect.boxx(boxx(c(-1,1),c(0,2),c(-1,1)),
                        boxx(c(0,3),c(0,1),c(0,4)))
    stopifnot(identical(C, boxx(c(0,1),c(0,1),c(0,1))))
    ## Should be empty (NULL)
    D <- intersect.boxx(boxx(c(-1,1),c(0,2),c(-1,1)),
                        boxx(c(0,3),c(0,1),c(0,4)), NULL)
    stopifnot(is.null(D))

    ## Tests for [.boxx with clip:
    ## Check ppp and ppx subset with clip are the same
    X <- cells
    WX <- shift(domain(X), vec = c(.5,.5))
    X2 <- X[WX, clip=TRUE]
    Y <- ppx(coords(X), domain = boxx(c(0,1),c(0,1)))
    WY <- shift(domain(Y), vec = c(.5,.5))
    Y2 <- Y[WY, clip=TRUE]
    stopifnot(all.equal(coords(X2), coords(Y2), check.attributes = FALSE))
    stopifnot(all.equal(domain(X2), as.owin(domain(Y2))))
  }

})
