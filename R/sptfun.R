#'
#'   sptfun.R
#'
#'   Support function of a (convex) set
#'   and related operations
#'
#'   $Revision: 1.6 $  $Date: 2026/02/09 08:00:55 $

SupportFun <- function(W, origin=c(0,0)) {
  stopifnot(is.owin(W))
  if(!missing(origin)) {
    origin <- interpretAsOrigin(origin, W)
    W <- shift(W, -origin)
  }
  V <- vertices(convexhull(W))
  f <- function(theta) { supportCalc(V, theta, "max") }
  return(f)
}

polarbody <- function(W, origin=c(0,0), npoly=512) {
  stopifnot(is.owin(W))
  if(!missing(origin)) {
    origin <- interpretAsOrigin(origin, W)
    W <- shift(W, -origin)
  }
  f <- SupportFun(W)
  theta <- seq(0, 2*pi, length=npoly)
  ra <- 1/f(theta)
  if(any(is.infinite(ra)))
    stop("Polar body has infinite radius; use a different origin point",
         call.=FALSE)
  PB <- owin(poly=list(x=ra * cos(theta), y=ra * sin(theta)),
            unitname=unitname(W))
  if(!missing(origin))
    PB <- shift(PB, origin)
  return(PB)
}

voronoiFlower <- function(W, origin=c(0,0), npoly=512) {
  stopifnot(is.owin(W))
  if(!missing(origin)) {
    origin <- interpretAsOrigin(origin, W)
    W <- shift(W, -origin)
  }
  f <- SupportFun(W)
  theta <- seq(0, 2*pi, length=npoly)
  ra <- 2 * f(theta)
  VF <- owin(poly=list(x=ra * cos(theta), y=ra * sin(theta)),
            unitname=unitname(W))
  if(!missing(origin))
    VF <- shift(VF, origin)
  return(VF)
}
  
FeretDiamFun <- function(W) {
  stopifnot(is.owin(W))
  V <- vertices(convexhull(W))
  f <- function(theta) { supportCalc(V, theta, "diffrange") }
  return(f)
}

minFeretDiam <- function(W) {
  stopifnot(is.owin(W))
  V <- vertices(convexhull(W))
  theta <- 2 * pi * seq(0, 1, length=2014)
  z <- supportCalc(V, theta, "diffrange")
  return(min(z))
}

maxFeretDiam <- function(W) {
  stopifnot(is.owin(W))
  V <- vertices(convexhull(W))
  theta <- 2 * pi * seq(0, 1, length=2014)
  z <- supportCalc(V, theta, "diffrange")
  return(max(z))
}

FeretBox <- function(W, theta=NULL) {
  stopifnot(is.owin(W))
  V <- vertices(convexhull(W))
  if(!is.null(theta)) {
    theta <- as.numeric(theta)
    theta <- c(theta, theta + pi/2) %% (2*pi)
  } else {
    theta <- 2 * pi * seq(0, 1, length=2048)
  }
  z <- supportCalc(V, theta, "diffrange")
  ## narrowest
  i0 <- which.min(z)
  theta0 <- theta[i0]
  ## 
  zhi.theta0 <- supportCalc(V, theta0, "max")
  zlo.theta0 <- supportCalc(V, theta0, "min")
  ## perpendicular to narrowest
  theta1 <- (theta0 + pi/2) %% (2 * pi)
  zhi.theta1 <- supportCalc(V, theta1, "max")
  zlo.theta1 <- supportCalc(V, theta1, "min")
  ## make box
  B <- owin(c(zlo.theta0, zhi.theta0), c(zlo.theta1, zhi.theta1))
  B <- rotate(B, theta0)
  unitname(B) <- unitname(W)
  return(B)
}

supportCalc <- function(V, theta, op=c("max", "min", "range", "diffrange")) {
  #' Calculate support function or similar
  #' Splitting large matrices into pieces
  op <- match.arg(op)

  x <- as.numeric(V$x)
  y <- as.numeric(V$y)
  nx <- length(x)

  dimt <- dim(theta)
  theta <- as.numeric(theta)
  nt <- length(theta)

  co <- cos(theta)
  si <- sin(theta)

  nm <- spatstat.options("maxmatrix")
  batchsize <- floor(nm/nt)
  if(nx <= batchsize) {
    M <- outer(x, co, "*") + outer(y, si, "*")
    z <- switch(op,
                max = apply(M, 2L, max),
                min = apply(M, 2L, min),
                range = apply(M, 2L, range),
                diffrange = apply(apply(M, 2L, range), 2L, diff)
                )
  } else {
    nbatches <- ceiling(nx/batchsize)
    for(ibatch in 1:nbatches) {
      batchstart <- 1L + (ibatch - 1L) * batchsize
      batchend   <- min(nx, ibatch * batchsize)
      j <- batchstart:batchend
      Mi <- outer(x[j], co, "*") + outer(y[j], si, "*")
      zi <- switch(op,
                   max = apply(Mi, 2L, max),
                   min = apply(Mi, 2L, min),
                   range = apply(Mi, 2L, range),
                   diffrange = apply(Mi, 2L, range) # postpone diff
                   )
      if(ibatch == 1) {
        z <- zi
      } else {
        z <- switch(op,
                    max = pmax(z, zi),
                    min = pmin(z, zi),
                    range = ,
                    diffrange = apply(rbind(z, zi), 2L, range))
      }
    }
    if(op == "diffrange") 
      z <- apply(z, 2L, diff)
  }
  if(!is.null(dimt)) 
    dim(z) <- c(dimt, nrow(z))
  return(z)
}


