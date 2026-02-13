#'
#'       minkowski.R
#' 
#'  Minkowski Sum and related operations
#'
#'  $Revision: 1.21 $ $Date: 2026/02/13 05:13:35 $


"%(+)%" <- MinkowskiSum <- function(A, B) { MinkowskiSumInternal(A,B) }

MinkowskiSumInternal <- local({

  MinkowskiSumInternal <- function(A, B, ..., do.specialcases=TRUE, p = NULL) {
    if(is.ppp(A)) {
      ## point pattern + something
      result <- UnionOfShifts(B, A)
    } else if(is.ppp(B)) {
      ## something + point pattern
      result <- UnionOfShifts(A, B)
    } else if(is.psp(A) && is.psp(B)) {
      ## segments + segments
      result <- UnionOfParallelograms(A,B)
    } else if(do.specialcases && is.rectangle(A) && is.rectangle(B)) {
      ## rectangle + rectangle
      result <- SumOfBoxes(A,B)
    } else {
      ## gather information
      simplyA <- simplyB <- convexA <- convexB <- FALSE
      if(owinA <- !is.psp(A)) {
        A <- as.polygonal(A)
        simplyA <- (length(A$bdry) == 1) # simply connected
        convexA <- is.convex(A)
      }
      if(owinB <- !is.psp(B)) {
        B <- as.polygonal(B)
        simplyB <- (length(B$bdry) == 1) # simply connected
        convexB <- is.convex(B)
      }
      ## convert each object to a list of polygons
      elementaryAA <- elementaryBB <- FALSE
      if(is.psp(A)) {
        ## convert line segments to flat polygons with 2 edges
        AA <- psp2poly(A)
        elementaryAA <- TRUE
      } else {
        ## convert window to list of polygons
        if(convexA || (do.specialcases && simplyA && convexB)) {
          ## A is a single polygon and we can use it as is
          AA <- list(A)
          elementaryAA <- FALSE
        } else {
          ## decompose A into triangles
          AA <- tiles(triangulate.owin(A))
          elementaryAA <- TRUE
        }
        AA <- lapply(AA, owin2singlepoly)
      }
      if(is.psp(B)) {
        ## convert line segments to flat polygons with 2 edges
        BB <- psp2poly(B)
        elementaryBB <- TRUE
      } else {
        ## convert window to list of polygons
        if(convexB || (do.specialcases && simplyB && convexA)) {
          ## B is a single polygon and we can use it as is
          BB <- list(B)
          elementaryBB <- FALSE
        } else {
          ## decompose B into triangles
          BB <- tiles(triangulate.owin(B))
          elementaryBB <- TRUE
        }
        BB <- lapply(BB, owin2singlepoly)
      }
      ## determine common resolution for polyclip operations
      FRAME <- SumOfBoxes(A, B)
      pdefault <- list(x0 = FRAME$xrange[1L],
                       y0 = FRAME$yrange[1L],
                       eps = mean(c(sidelengths(Frame(A)),
                                    sidelengths(Frame(B))))/2^30)
      p <- resolve.defaults(p, pdefault)
      ## compute Minkowski sums of pairs of CONVEX pieces
      result <- NULL
      if(do.specialcases && elementaryAA && elementaryBB) {
        ## each piece in AA and in BB is a triangle or line segment; evaluate directly
        for(a in AA) {
          partial.a <- NULL
          for(b in BB) {
            x <- as.numeric(outer(a$x, b$x, "+"))
            y <- as.numeric(outer(a$y, b$y, "+"))
            h <- rev(chull(x,y))
            if(length(h) >= 3) {
              contrib.ab <- owin(poly=list(x=x[h], y=y[h]), check=FALSE) 
              partial.a <- union.owin(partial.a, contrib.ab, p=p)
            }
          }
          if(!is.null(partial.a)) 
            result <- union.owin(result, partial.a, p=p)
        }
      } else {
        ## general case - use polyclip
        x0 <- p$x0
        y0 <- p$y0
        eps <- p$eps
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
      if(!is.rectangle(result))
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

  owin2singlepoly <- function(w) { w$bdry[[1]] }
  
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
    nA <- nrow(eA)
    nB <- nrow(eB)
    ## determine common resolution for polyclip operations
    eps <- mean(c(sidelengths(Frame(A)), sidelengths(Frame(B))))/2^30
    p <- list(eps=eps)
    ## form the Minkowski sum of each pair of segments
    flat <- list()
    interior <- NULL
    for(i in seq_len(nA)) {
      for(j in seq_len(nB)) {
        x <- as.numeric(outer(eA[i, c(1,3)], eB[j, c(1,3)], "+"))
        y <- as.numeric(outer(eA[i, c(2,4)], eB[j, c(2,4)], "+"))
        h <- rev(chull(x,y))
        xy <- list(x=x[h], y=y[h])
        if(length(h) >= 3) {
          ## polygon with interior
          g <- owin(poly=xy, check=FALSE)
          interior <- union.owin(interior, g, p=p)
        } else if(length(h) == 2) {
          ## line segment: store it
          flat <- append(flat, list(xy))
        } 
      }
    }
    if(!is.null(interior)) {
      ## fill holes due to numerical error 
      amin <- sqrt(.Machine$double.eps) * area(interior)
      result <- fillholes.owin(interior, amin)
    } else {
      ## all segments were parallel
      ## result is a psp
      xx <- unlist(lapply(flat, getelement, name="x"))
      yy <- unlist(lapply(flat, getelement, name="y"))
      ends <- cbind(xx[c(TRUE,  FALSE)],
                    yy[c(TRUE,  FALSE)],
                    xx[c(FALSE, TRUE)],
                    yy[c(FALSE, TRUE)])
      result <- as.psp(ends, window=SumOfBoxes(A,B), check=FALSE)
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
  
  MinkowskiSumInternal
})

dilationAny <- function(A, B) { MinkowskiSumInternal(A, reflect(B)) }

"%(-)%" <- erosionAny <- function(A, B) {
  D <- Frame(A)
  Dplus <- grow.rectangle(D, 0.1 * shortside(D))
  Ac <- complement.owin(A, Dplus)
  AcB <- MinkowskiSumInternal(Ac, reflect(B))
  if(is.subset.owin(D, AcB))
    return(emptywindow(D))
  C <- complement.owin(AcB[Dplus], Dplus)[D]
  return(C)
}
