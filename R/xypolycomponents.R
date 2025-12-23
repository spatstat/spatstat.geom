#'
#'  xypolycomponents.R
#'
#'  Relations between component polygons of a window
#'
#'  $Revision: 1.5 $ $Date: 2025/12/23 01:51:20 $

xypolycomponents <- function(W) {
  stopifnot(is.polygonal(W))
  A <- xypolyNesting(W)
  B <- W$bdry
  uW <- unitname(W)
  solapply(A$components, xypolysub2owin, B=B, unitname=uW)
}

xypoly2owin <- function(p, check=FALSE, ...) {
  rescue.rectangle(owin(poly=p, check=check, ...))
}

xypolysub2owin <- function(i, B, ...) { xypoly2owin(B[i], ...) }

xypolyNesting <- function(W) {
  B <- W$bdry
  nB <- length(B)
  M <- matrix(FALSE, nB, nB)
  diag(M) <- FALSE
  h <- sapply(B, is.hole.xypolygon)
  if(nB > 1 && any(h)) {
    G <- B
    G[h] <- lapply(B[h], reverse.xypolygon)
    G <- lapply(G, xypoly2owin)
    for(i in 2:nB) {
      Gi <- G[[i]]
      for(j in 1:(i-1)) {
        Gj <- G[[j]]
        M[i,j] <- Mij <- is.subset.owin(Gi, Gj)
        ## if M[i,j] = TRUE then M[j,i] = FALSE so take no action 
        if(!Mij) 
          M[j,i] <- is.subset.owin(Gj, Gi)
      }
    }
  }
  ## Count[i] is the number of contours that enclose B[[i]]
  Count <- rowSums(M)
  p <- !h
  ## count is odd iff contour is hole
  if(!all((Count %% 2 == 0) != h))
    stop("Contour counts do not agree with hole status", call.=FALSE)
  ## Order[i] is the number of positive contours that enclose B[[i]]
  ## including B[[i]] itself. 
  Order <- p + rowSums(M[, p, drop=FALSE])
  ## identify connected components
  Compo <- xypolycoco(M, h)
  result <- list(hole=h, subset=M, count=Count, order=Order, components=Compo)
  return(result)
}

## identify connected components 

xypolycoco <- function(M, h, id) {
  Count <- rowSums(M)
  p <- !h
  Order <- p + rowSums(M[, p, drop=FALSE])  
  if(missing(id)) id <- seq_along(h)
  ## identify exterior contours
  exterior <- (Count == 0)
  nExt <- sum(exterior)
  if(nExt == 0)
    stop("No exterior contour found", call.=FALSE)
  if(nExt > 1) {
    ## several exterior contours
    ## (1) sanity check:
    if(any(rowSums(M[!exterior, exterior, drop=FALSE]) > 1))
      stop("Subset relations are paradoxical", call.=FALSE)
    ## (2) segregate interior contours according to their exterior contour
    result <- list()
    for(iext in which(exterior)) {
      jsub <- c(iext, which(!exterior & rowSums(M[, iext, drop=FALSE]) > 0))
      Msub <- M[jsub, jsub, drop=FALSE]
      hsub <- h[jsub]
      idsub <- id[jsub]
      result <- c(result, xypolycoco(Msub, hsub, idsub))
    }
    return(result)
  }
  ## single exterior contour: descend
  outerlayer <- (Count <= 1)
  result <- list(id[outerlayer])
  if(!all(outerlayer)) {
    interior <- !outerlayer
    ## strip off outer layer and recurse
    Msub <- M[interior, interior, drop=FALSE]
    hsub <- h[interior]
    idsub <- id[interior]
    result <- c(result, xypolycoco(Msub, hsub, idsub))
  }
  return(result)
}

