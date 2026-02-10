#'
#'       minkowski.R
#' 
#'  Minkowski Sum and related operations
#'
#'  $Revision: 1.15 $ $Date: 2026/02/10 03:58:59 $


"%(+)%" <- MinkowskiSum <- local({

  MinkowskiSum <- function(A, B) {
    if(is.ppp(A)) {
      ## point pattern + something
      result <- UnionOfShifts(B, A)
    } else if(is.ppp(B)) {
      ## something + point pattern
      result <- UnionOfShifts(A, B)
    } else if(is.psp(A) && is.psp(B)) {
      ## segments + segments
      result <- UnionOfParallelograms(A,B)
    } else if(is.rectangle(A) && is.rectangle(B)) {
      result <- SumOfBoxes(A,B)
    } else {
      ## convert each object to a list of CONVEX polygons
      if(is.psp(A)) {
        AA <- psp2poly(A)
        simpleA <- TRUE
      } else {
        A <- as.polygonal(A)
        if(is.convex(A)) {
          AA <- list(A)
          simpleA <- FALSE
        } else {
          AA <- tiles(triangulate.owin(A))
          simpleA <- TRUE
        }
        AA <- lapply(AA, owin2poly)
      }
      if(is.psp(B)) {
        BB <- psp2poly(B)
        simpleB <- TRUE
      } else {
        B <- as.polygonal(B)
        if(is.convex(B)) {
          BB <- list(B)
          simpleB <- FALSE
        } else {
          BB <- tiles(triangulate.owin(B))
          simpleB <- TRUE
        }
        BB <- lapply(BB, owin2poly)
      }
      ## determine common resolution for polyclip operations
      FRAME <- SumOfBoxes(A, B)
      x0 <- FRAME$xrange[1L]
      y0 <- FRAME$yrange[1L]
      eps <- mean(c(sidelengths(Frame(A)), sidelengths(Frame(B))))/2^30
      p <- list(eps=eps, x0=x0, y0=y0)
      ## compute Minkowski sums of pairs of CONVEX pieces
      result <- NULL
      if(simpleA && simpleB) {
        ## each piece is a triangle or line segment; evaluate directly
        for(a in AA) {
          partial.a <- NULL
          for(b in BB) {
            x <- as.numeric(outer(a$x, b$x, "+"))
            y <- as.numeric(outer(a$y, b$y, "+"))
            h <- rev(chull(x,y))
            contrib.ab <- owin(poly=list(x=x[h], y=y[h]), check=FALSE) 
            partial.a <- union.owin(partial.a, contrib.ab, p=p)
          }
          result <- union.owin(result, partial.a, p=p)
        }
      } else {
        ## general case - use polyclip
        for(a in AA) {
          partial.a <- NULL
          for(b in BB) {
            contrib.ab <- polyclip::polyminkowski(a, b, x0=x0, y0=y0, eps=eps)
            contrib.ba <- polyclip::polyminkowski(b, a, x0=x0, y0=y0, eps=eps)
            ## these must be convex
            contrib.ab <- fillholes.owin(poly2owin(contrib.ab), Inf)
            contrib.ba <- fillholes.owin(poly2owin(contrib.ba), Inf)
            ## union
            partial.a <- union.owin(partial.a, contrib.ab, contrib.ba, p=p)
          }
          result <- union.owin(result, partial.a, p=p)
        }
      }
      result <- rescue.rectangle(result)
      Frame(result) <- FRAME
    }
    ## resolve unitname
    un <- list(unitname(A), unitname(B))
    un <- unique(un[!sapply(un, is.vanilla)])
    if(length(un) == 1)
      unitname(result) <- un[[1L]]
    return(result)
  }

  ## processing helper functions

  poly2owin <- function(z) owin(poly=z, check=FALSE)

  owin2poly <- function(w) { w$bdry }
  
  psp2poly <- function(X) apply(as.matrix(X$ends), 1, seg2poly)

  seg2poly <- function(z) with(as.list(z), list(x=c(x0, x1, x0), y=c(y0,y1,y0)))

  ##
  UnionOfShifts <- function(X, V) {
    #' compute the union or superposition of copies of X by vectors in V
    v <- as.matrix(coords(V))
    n <- nrow(v)
    Y <- vector(mode="list", length=n)
    for(i in seq_len(n)) 
      Y[[i]] <- shift(X, v[i,])
    Y <- as.solist(Y)
    if(is.owin(X)) {
      Z <- union.owin(Y)
    } else {
      #' X is a pattern of objects in a window
      W <- MinkowskiSum(Window(X), Window(V))
      Z <- superimpose(Y, W=W)
    }
    return(Z)
  }

  UnionOfParallelograms <- function(A,B) {
    ## A and B are segment patterns
    eA <- as.matrix(A$ends)
    eB <- as.matrix(B$ends)
    ## determine common resolution for polyclip operations
    eps <- mean(c(sidelengths(Frame(A)), sidelengths(Frame(B))))/2^30
    p <- list(eps=eps)
    ## form the Minkowski sum of each pair of segments
    result <- NULL
    for(i in seq_len(nrow(eA))) {
      for(j in seq_len(nrow(eB))) {
        x <- as.numeric(outer(eA[i, c(1,3)], eB[j, c(1,3)], "+"))
        y <- as.numeric(outer(eA[i, c(2,4)], eB[j, c(2,4)], "+"))
        h <- rev(chull(x,y))
        g <- owin(poly=list(x=x[h], y=y[h]), check=FALSE) 
        result <- union.owin(result, g, p=p)
      }
    }
    return(result)
  }

  SumOfBoxes <- function(A, B) {
    A <- Frame(A)
    B <- Frame(B)
    xr <- A$xrange + B$xrange
    yr <- A$yrange + B$yrange
    result <- owin(xr, yr)
    return(result)
  }
  
  MinkowskiSum
})

dilationAny <- function(A, B) { MinkowskiSum(A, reflect(B)) }

"%(-)%" <- erosionAny <- function(A, B) {
  D <- Frame(A)
  Dplus <- grow.rectangle(D, 0.1 * shortside(D))
  Ac <- complement.owin(A, Dplus)
  AcB <- MinkowskiSum(Ac, reflect(B))
  if(is.subset.owin(D, AcB))
    return(emptywindow(D))
  C <- complement.owin(AcB[Dplus], Dplus)[D]
  return(C)
}
