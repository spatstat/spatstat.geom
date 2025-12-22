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
##
##  tests/segments.R
##   Tests of psp class and related code
##                      [SEE ALSO: tests/xysegment.R]
##
##  $Revision: 1.33 $  $Date: 2022/05/22 08:39:47 $


local({

  if(ALWAYS) { # depends on platform
    ## pointed out by Jeff Laake
    W <- owin()
    X <- psp(x0=.25,x1=.25,y0=0,y1=1,window=W)
    X[W]
  }

  X <- psp(runif(10),runif(10),runif(10),runif(10), window=owin())

  if(FULLTEST) {
    Z <- as.mask.psp(X)
    Z <- pixellate(X)
  }

  if(ALWAYS) { # platform dependent
    ## add short segment
    Shorty <- psp(0.5, 0.6, 0.5001, 0.6001, window=Window(X))
    XX <- superimpose(X[1:5], Shorty, X[6:10])
    ZZ <- as.mask.psp(XX)
    ZZ <- pixellate(XX)
  }

  if(FULLTEST) {
    #' misc
    PX <- periodify(X, 2)
  }

  if(ALWAYS) { # C code
    ## tests of pixellate.psp -> seg2pixL
    ns <- 50
    out <- numeric(ns)
    for(i in 1:ns) {
      X <- psp(runif(1), runif(1), runif(1), runif(1), window=owin())
      len <- lengths_psp(X)
      dlen <- sum(pixellate(X)$v)
      out[i] <- if(len > 1e-7) dlen/len else 1
    }
    if(diff(range(out)) > 0.01) stop(paste(
                                  "pixellate.psp test 1: relative error [",
                                  paste(diff(range(out)), collapse=", "),
                                  "]"))

    ## Michael Sumner's test examples
    set.seed(33)
    n <- 2001
    co <- cbind(runif(n), runif(n))
    ow <- owin()
    X <- psp(co[-n,1], co[-n,2], co[-1,1], co[-1,2], window=ow)
    s1 <- sum(pixellate(X))
    s2 <- sum(lengths_psp(X))
    if(abs(s1 - s2)/s2 > 0.01) {
      stop(paste("pixellate.psp test 2:",
                 "sum(pixellate(X)) = ", s1,
                 "!=", s2, "= sum(lengths_psp(X))"))
    }

    wts <- 1/(lengths_psp(X) * X$n)
    s1 <- sum(pixellate(X, weights=wts))
    if(abs(s1-1) > 0.01) {
      stop(paste("pixellate.psp test 3:",
                 "sum(pixellate(X, weights))=", s1,
                 " (should be 1)"))
    }
    
    X <- psp(0, 0, 0.01, 0.001, window=owin())
    s1 <- sum(pixellate(X))
    s2 <- sum(lengths_psp(X))
    if(abs(s1 - s2)/s2 > 0.01) {
      stop(paste("pixellate.psp test 4:",
                 "sum(pixellate(X)) = ", s1,
                 "!=", s2, "= sum(lengths_psp(X))"))
    }

    X <- psp(0, 0, 0.001, 0.001, window=owin())
    s1 <- sum(pixellate(X))
    s2 <- sum(lengths_psp(X))
    if(abs(s1 - s2)/s2 > 0.01) {
      stop(paste("pixellate.psp test 5:",
                 "sum(pixellate(X)) = ", s1,
                 "!=", s2, "= sum(lengths_psp(X))"))
    }
  }

  if(FULLTEST) {
    #' cases of superimpose.psp
    A <- as.psp(matrix(runif(40), 10, 4), window=owin())
    B <- as.psp(matrix(runif(40), 10, 4), window=owin())
    superimpose(A, B, W=ripras)
    superimpose(A, B, W="convex")
  }

  if(FULLTEST) {
    #' as.psp.data.frame
    df <- as.data.frame(matrix(runif(40), ncol=4))
    A <- as.psp(df, window=square(1))
    colnames(df) <- c("x0","y0","x1","y1")
    df <- cbind(df, data.frame(marks=1:nrow(df)))
    B <- as.psp(df, window=square(1))
    colnames(df) <- c("xmid", "ymid", "length", "angle", "marks")
    E <- as.psp(df, window=square(c(-1,2)))
    G <- E %mark% factor(sample(letters[1:3], nsegments(E), replace=TRUE))
    H <- E %mark% runif(nsegments(E))
  
    #' print and summary methods
    A
    B
    E
    G
    H
    summary(B)
    summary(G)
    summary(H)
    M <- B
    marks(M) <- data.frame(id=marks(B), len=lengths_psp(B))
    M
    summary(M)
    subset(M, select=len)

    #' plot method cases  
    spatstat.options(monochrome=TRUE)
    plot(B)
    plot(G)
    plot(M)
    spatstat.options(monochrome=FALSE)
    plot(B)
    plot(G)
    plot(M)
    #' misuse of 'col' argument - several cases
    plot(G, col="grey") # discrete
    plot(B, col="grey") 
    plot(unmark(B), col="grey") 
    plot(M, col="grey") 
  
    #' miscellaneous class support cases
    marks(M) <- marks(M)[1,,drop=FALSE]
    
    #' undocumented  
    as.ppp(B)
  }

  if(ALWAYS) { # C code
    #' segment crossing code
    X <- psp(runif(30),runif(30),runif(30),runif(30), window=owin())
    A <- selfcut.psp(X, eps=1e-11)
    B <- selfcut.psp(X[1])
    #' 
    Y <- psp(runif(30),runif(30),runif(30),runif(30), window=owin())
    Z <- edges(letterR)[c(FALSE,TRUE)]
    spatstat.options(selfcrossing.psp.useCall=FALSE, crossing.psp.useCall=FALSE)
    A <- selfcrossing.psp(X)
    B <- selfcrossing.psp(Z)
    D <- crossing.psp(X,Y,details=TRUE)
    spatstat.options(selfcrossing.psp.useCall=TRUE, crossing.psp.useCall=TRUE)
    A <- selfcrossing.psp(X)
    B <- selfcrossing.psp(Z)
    D <- crossing.psp(X,Y,details=TRUE)
    reset.spatstat.options()
  }

  if(FULLTEST) {
    #' geometry
    m <- data.frame(A=1:10, B=letters[1:10])
    X <- psp(runif(10), runif(10), runif(10), runif(10), window=owin(), marks=m)
    Z <- rotate(X, angle=pi/3, centre=c(0.5, 0.5))
    Y <- endpoints.psp(X, which="lower")
    Y <- endpoints.psp(X, which="upper")
    Y <- endpoints.psp(X, which="right")
    U <- flipxy(X)
  }
  
  if(ALWAYS) {
    ## nnfun.psp
    P <- psp(runif(10), runif(10), runif(10), runif(10),
             window=square(1), marks=runif(10))
    f <- nnfun(P)
    f <- nnfun(P, value="mark")
    d <- domain(f)
    Z <- as.im(f)
  }

})

reset.spatstat.options()





#'
#'   spatstat.geom/tests/sidefx.R
#'
#' Test whether plot(do.plot=FALSE) has no side effects on graphics system
#'
#'  $Revision: 1.2 $  $Date: 2025/12/22 02:25:23 $

local({
  if(FULLTEST) {
    ## test whether a graphics device has been started
    deviceExists <- function() { length(dev.list()) != 0 }
    ## check whether executing 'expr' causes creation of a graphics device
    chk <- function(expr) {
      ename <- sQuote(deparse(substitute(expr)))
      if(deviceExists()) {
        ## try switching off the graphics
        graphics.off()
        if(deviceExists()) {
          ## Dang
          warning(paste("Cannot check", ename, 
                        "as a graphics device already exists"),
                  call.=FALSE)
          return(FALSE)
        }
      }
      eval(expr)
      if(deviceExists()) {
        stop(paste("Evaluating", ename, 
                   "caused a graphics device to be started"),
             call.=FALSE)
      }
      return(TRUE)
    }
    
    ## windows
    chk(plot(owin(), do.plot=FALSE))
    chk(plot(letterR, do.plot=FALSE))
    chk(plot(as.mask(letterR), do.plot=FALSE))
    ## point patterns
    chk(plot(cells, do.plot=FALSE))
    chk(plot(chorley, do.plot=FALSE))
    chk(plot(chorley, legend=FALSE, do.plot=FALSE))
    ## line segment patterns
    chk(plot(edges(letterR), do.plot=FALSE))
    ## images
    chk(plot(bei.extra$elev, do.plot=FALSE))
    chk(plot(bei.extra$elev, ribbon=FALSE, do.plot=FALSE))
    ## solist/imlist
    chk(plot(bei.extra, do.plot=FALSE))
    chk(plot(solist(cells, redwood), do.plot=FALSE))
  }
})
#'
#'     tests/simplepan.R
#'
#'   Tests of user interaction in simplepanel
#'   Handled by spatstatLocator()
#'
#'   $Revision: 1.3 $  $Date: 2020/05/01 09:59:59 $
#'

if(ALWAYS) {  # may depend on platform
local({
  ## Adapted from example(simplepanel)
  ## make boxes
  outerbox <- owin(c(0,4), c(0,1))
  buttonboxes <- layout.boxes(outerbox, 4, horizontal=TRUE, aspect=1)
  ## make environment containing an integer count
  myenv <- new.env()
  assign("answer", 0, envir=myenv)
  ## what to do when finished: return the count.
  myexit <- function(e) { return(get("answer", envir=e)) }
  ## button clicks
  ## decrement the count
  Cminus <- function(e, xy) {
    ans <- get("answer", envir=e)
    assign("answer", ans - 1, envir=e)
    return(TRUE)
  }
  ## display the count (clicking does nothing)
  Cvalue <- function(...) { TRUE }
  ## increment the count
  Cplus <- function(e, xy) {
    ans <- get("answer", envir=e)
    assign("answer", ans + 1, envir=e)
    return(TRUE)
  }
  ## 'Clear' button
  Cclear <- function(e, xy) {
    assign("answer", 0, envir=e)
    return(TRUE)
  }
  ## quit button
  Cdone <- function(e, xy) { return(FALSE) }

  myclicks <- list("-"=Cminus,
                   value=Cvalue,
                   "+"=Cplus,
                   done=Cdone)
  ## redraw the button that displays the current value of the count
  Rvalue <- function(button, nam, e) {
     plot(button, add=TRUE)
     ans <- get("answer", envir=e)
     text(centroid.owin(button), labels=ans)
     return(TRUE)
  }
  ## make the panel
  P <- simplepanel("Counter",
                   B=outerbox, boxes=buttonboxes,
                   clicks=myclicks,
                   redraws = list(NULL, Rvalue, NULL, NULL),
                   exit=myexit, env=myenv)
  ## queue up a sequence of inputs
  boxcentres <- do.call(concatxy, unname(lapply(buttonboxes[c(3,3,1,3,2,4)],
                                                centroid.owin)))
  spatstat.utils::queueSpatstatLocator(boxcentres$x, boxcentres$y)
  ## go
  run.simplepanel(P)
})
}
#
#  tests/splitpea.R
#
#  Check behaviour of split.ppp etc
#
#  Thanks to Marcelino de la Cruz
#
#  $Revision: 1.17 $  $Date: 2021/04/15 06:19:51 $
#

local({
  W <- square(8)
  X <- ppp(c(2.98, 4.58, 7.27, 1.61, 7.19),
           c(7.56, 5.29, 5.03, 0.49, 1.65),
           window=W, check=FALSE)
  Z <- quadrats(W, 4, 4)
  Yall <- split(X, Z, drop=FALSE)
  Ydrop <- split(X, Z, drop=TRUE)

  if(ALWAYS) { # may depend on platform
    P <- Yall[[1]]
    if(!all(inside.owin(P$x, P$y, P$window)))
      stop("Black hole detected when drop=FALSE")
    P <- Ydrop[[1]]
    if(!all(inside.owin(P$x, P$y, P$window)))
      stop("Black hole detected when drop=TRUE")
    Ydrop[[1]] <- P[1]
    split(X, Z, drop=TRUE) <- Ydrop
  }

  ## test NA handling
  Zbad <- quadrats(square(4), 2, 2)
  Ybdrop <- split(X, Zbad, drop=TRUE)
  Yball  <- split(X, Zbad, drop=FALSE)

  if(FULLTEST) {
    ## other bugs/ code blocks in split.ppp, split<-.ppp, [<-.splitppp
    flog <- rep(c(TRUE,FALSE), 21)
    fimg <- as.im(dirichlet(runifrect(5, Window(cells))), dimyx=32)
    A <- split(cells, flog)
    B <- split(cells, square(0.5))
    D <- split(cells, fimg)
    E <- split(cells, logical(42), drop=TRUE)
    Cellules <- cells
    split(Cellules, flog) <- solapply(A, rjitter)
    split(Cellules, fimg) <- solapply(D, rjitter)
    D[[2]] <- rjitter(D[[2]])
    Funpines <- finpines
    marks(Funpines)[,"diameter"] <- factor(marks(Funpines)[,"diameter"])
    G <- split(Funpines)
    H <- split(Funpines, "diameter")
    split(Funpines) <- solapply(G, rjitter)
    split(Funpines, "diameter") <- solapply(H, rjitter)

    ## From Marcelino
    set.seed(1)
    W<- square(10) # the big window
    ## puntos<- rpoispp(0.5, win=W)
    puntos<- runifrect(rpois(1, 0.5 * area(W)), win=W)
    r00 <- letterR
    r05 <- shift(letterR,c(0,5))
    r50 <- shift(letterR,c(5,0))
    r55 <- shift(letterR,c(5,5))
    tessr4 <- tess(tiles=list(r00, r05,r50,r55))
    puntosr4 <- split(puntos, tessr4, drop=TRUE)
    split(puntos, tessr4, drop=TRUE) <- puntosr4

    ## More headaches with mark format
    A <- runifrect(10)
    B <- runifrect(10)
    AB <- split(superimpose(A=A, B=B))

    #' check that split<- respects ordering where possible
    X <- amacrine
    Y <- split(X)
    split(X) <- Y
    stopifnot(identical(X, amacrine))

    #' split.ppx
    df <- data.frame(x=runif(4),y=runif(4),t=runif(4),
                     age=rep(c("old", "new"), 2),
                     mineral=factor(rep(c("Au","Cu"), each=2),
                                    levels=c("Au", "Cu", "Pb")),
                     size=runif(4))
    X <- ppx(data=df, coord.type=c("s","s","t","m", "m","m"))
    Y <- split(X, "age")
    Y <- split(X, "mineral", drop=TRUE)
    Y <- split(X, "mineral")
    print(Y)
    print(summary(Y))
    Y[c(TRUE,FALSE,TRUE)]
    Y[1:2]
    Y[3] <- Y[1]
  }
})

##
## tests/symbolmaps.R
##
##   Quirks associated with symbolmaps, etc.
##
## $Revision: 1.6 $ $Date: 2024/04/08 04:22:25 $

if(FULLTEST) {
local({
  set.seed(100)
  X <- runifrect(8)

  ## symbolmap for numeric values
  g1 <- symbolmap(range=c(0,100), size=function(x) x/50)
  invoke.symbolmap(g1, 50, x=numeric(0), y=numeric(0), add=TRUE)
  plot(g1, labelmap=100)
  ## symbolmap for discrete categories
  g2 <- symbolmap(inputs=letters[1:5], chars=1:5)
  invoke.symbolmap(g2, "a", x=numeric(0), y=numeric(0), add=TRUE)
  plot(g2)
  ## constant/trivial
  a <- symbolmap(pch=16)
  print(a)
  plot(a)
  symbolmapdomain(a)
  b <- symbolmap()
  print(b)
  ## graphical arguments with mixed types (function, constant)
  f <- function(x) { ifelse(x %in% letters[1:3], "circles", "squares")}
  g3 <- symbolmap(inputs=letters[1:5], size=0.7, shape=f)
  invoke.symbolmap(g3, "a", x=numeric(0), y=numeric(0), add=TRUE)
  plot(g3)

  ## textureplot
  V <- as.im(dirichlet(X))
  tmap <- textureplot(V)
  textureplot(V, textures=tmap, legend=TRUE, leg.side="left")
  textureplot(V, leg.side="bottom")
  textureplot(V, leg.side="top")
  ## spacing too large for tiles - upsets various pieces of code
  textureplot(V, spacing=2)

  ## plot.texturemap
  plot(tmap, vertical=TRUE)
  plot(tmap, vertical=TRUE, xlim=c(0,1))
  plot(tmap, vertical=TRUE, ylim=c(0,1))
  plot(tmap, vertical=FALSE, xlim=c(0,1))
  plot(tmap, vertical=FALSE, ylim=c(0,1))

  ## infrastructure
  plan.legend.layout(owin(), side="top", started=TRUE)
})
}
