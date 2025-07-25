#
# connected.R
#
# connected component transform
#
#    $Revision: 1.30 $  $Date: 2025/07/15 03:22:49 $
#
# Interpreted code for pixel images by Julian Burgos <jmburgos@u.washington.edu>
# Rewritten in C by Adrian Baddeley
#
# Code for point patterns by Adrian Baddeley

connected <- function(X, ...) {
  UseMethod("connected")
}

connected.im <- function(X, ..., background=NA, method="C", connect=8) {
  if(!is.na(background)) {
    W <- solutionset(X != background)
  } else if(X$type == "logical") {
    W <- solutionset(X)
  } else {
    warning("Assuming background = NA, foreground = other values", call.=FALSE)
    W <- as.owin(X)
  }
  connected.owin(W, method=method, ..., connect=connect)
}

connected.owin <- function(X, ..., polygonal=FALSE, method="C", connect=8) {
  if(polygonal) {
    P <- as.polygonal(X)
    A <- xypolycomponents(P)
    W <- if(is.mask(X)) P else X
    result <- tess(tiles=A, window=W)
    return(result)
  }
  method <- pickoption("algorithm choice", method,
                       c(C="C", interpreted="interpreted"))
  if(!missing(connect)) {
    check.1.integer(connect)
    if(!any(connect == c(4,8)))
      stop("'connect' should be 4 or 8")
  }
  # convert X to binary mask
  X <- as.mask(X, ...)
  #     
  Y <- X$m
  nr <- X$dim[1L]
  nc <- X$dim[2L]

  if(method == "C") {
################ COMPILED CODE #########################
# Pad border with FALSE
    M <- rbind(FALSE, Y, FALSE)
    M <- cbind(FALSE, M, FALSE)
    # assign unique label to each foreground pixel 
    L <- M
    L[M] <- seq_len(sum(M))
    L[!M] <- 0
    ## resolve labels
    if(typeof(L) == "double") {
      ## Labels are numeric (not integer)
      ## This can occur if raster is really large
      if(connect == 8) {
        ## 8-connected
        z <- .C(SG_coco8dbl,
                mat=as.double(t(L)),
                nr=as.integer(nr),
                nc=as.integer(nc),
                PACKAGE="spatstat.geom")
      } else {
        ## 4-connected
        z <- .C(SG_coco4dbl,
                mat=as.double(t(L)),
                nr=as.integer(nr),
                nc=as.integer(nc),
                PACKAGE="spatstat.geom")
      }
    } else {
      ## Labels are integer
      if(connect == 8) {
        ## 8-connected
        z <- .C(SG_coco8int,
                mat=as.integer(t(L)),
                nr=as.integer(nr),
                nc=as.integer(nc),
                PACKAGE="spatstat.geom")
      } else {
        ## 4-connected
        z <- .C(SG_coco4int,
                mat=as.integer(t(L)),
                nr=as.integer(nr),
                nc=as.integer(nc),
                PACKAGE="spatstat.geom")
      }
    }
    # unpack
    Z <- matrix(z$mat, nr+2, nc+2, byrow=TRUE)
  } else {
################ INTERPRETED CODE #########################
# by Julian Burgos
#  
# Pad border with zeros
    padY <- rbind(0, Y, 0)
    padY <- cbind(0, padY, 0)
    # Initialise 
    Z <- matrix(0, nrow(padY), ncol(padY))
    currentlab <- 1L
    todo <- as.vector(t(Y))
    equiv <- NULL
    ## extension by Adrian
    nprev <- as.integer(connect/2)
    
    # ........ main loop ..........................
    while(any(todo)){
      # pick first unresolved pixel
      one <- which(todo)[1L]
      onerow <- ceiling(one/nc)
      onecol <- one -((onerow-1L)*nc)
      parow=onerow+1L # Equivalent rows & column in padded matrix
      pacol=onecol+1L
      ## Examine four previously scanned neighbors
      ## (use padded matrix to avoid edge issues)
      nbrs <- if(connect == 8) {
                rbind(c(parow-1L,pacol-1L),
                      c(parow-1L,pacol),
                      c(parow,  pacol-1L),
                      c(parow-1L,pacol+1L))
              } else {
                rbind(c(parow-1L,pacol),
                      c(parow,  pacol-1L))
              }
      px <- sum(padY[nbrs])
      if (px==0){
        # no neighbours: new component
        Z[parow,pacol] <- currentlab
        currentlab <- currentlab+1L
        todo[one] <- FALSE
      } else if(px==1L) {
        # one neighbour: assign existing label
        labs <- unique(Z[nbrs], na.rm=TRUE)
        labs <- labs[labs != 0]
        Z[parow,pacol] <- labs[1L]
        currentlab <- max(Z)+1L
        todo[one] <- FALSE
      } else {
        # more than one neighbour: possible merger of labels
        labs <- unique(Z[nbrs], na.rm=TRUE)
        labs <- labs[labs != 0]
        labs <- sort(labs)
        equiv <- rbind(equiv,c(labs,rep.int(0,times=nprev-length(labs))))
        Z[parow,pacol] <- labs[1L]
        currentlab <- max(Z)+1L
        todo[one] <- FALSE
      }
    }
    # ........... end of loop ............
    # Resolve equivalences ................

    if(length(equiv)>1L){
      merges <- (equiv[,2L] > 1L)
      nmerge <- sum(merges)
      if(nmerge==1L)
        equiv <- equiv[which(merges), , drop=FALSE]
      else if(nmerge > 1L) {
        relevant <- (equiv[,2L] > 0)
        equiv <- equiv[relevant, , drop=FALSE]
        equiv <- equiv[fave.order(equiv[,1L]),]
      }
      for (i in 1:nrow(equiv)){
        current <- equiv[i, 1L]
        for (j in 2:nprev){
          twin <- equiv[i,j]
          if (twin>0){
            # Change labels matrix
            Z[which(Z==twin)] <- current
            # Update equivalence table
            equiv[which(equiv==twin)] <- current
          }
        }
      }
    }
  }

  ########### COMMON CODE ############################
    
  # Renumber labels sequentially
  mapped <- (Z != 0)
  usedlabs <- sortunique(as.vector(Z[mapped]))
  nlabs <- length(usedlabs)
  labtable <- cumsum(seq_len(max(usedlabs)) %in% usedlabs)
  Z[mapped] <- labtable[Z[mapped]]

  # banish zeroes
  Z[!mapped] <- NA
  
  # strip borders
  Z <- Z[2:(nrow(Z)-1L),2:(ncol(Z)-1L)]
  # dress up 
  Z <- im(factor(Z, levels=1:nlabs),
          xcol=X$xcol, yrow=X$yrow, unitname=unitname(X))
  return(Z)
}

connected.ppp <- connected.pp3 <- function(X, R, ...) {
  methodname <- if(is.ppp(X)) "connected.ppp" else
                if(is.pp3(X)) "connected.pp3" else 
                stopifnot(is.ppp(X) || is.pp3(X))
  check.1.real(R, paste("In", methodname))
  stopifnot(R >= 0)
  internal <- resolve.1.default("internal", list(...), list(internal=FALSE))
  nv <- npoints(X)
  cl <- closepairs(X, R, what="indices")
  lab0 <- cocoEngine(nv, cl$i - 1L, cl$j - 1L, methodname)
  if(internal)
    return(lab0)
  lab <- lab0 + 1L
  # Renumber labels sequentially 
  lab <- as.integer(factor(lab))
  # Convert labels to factor
  lab <- factor(lab)
  # Apply to points
  Y <- X %mark% lab
  return(Y)
}

cocoEngine <- function(nv, ie, je, algoname="connectedness algorithm") {
  #' no checks
  #' ie, je are 0-based indices (range between 0 and nv-1)
  ne <- length(ie)
  zz <- .C(SG_cocoGraph,
           nv=as.integer(nv),
           ne=as.integer(ne),
           ie=as.integer(ie),
           je=as.integer(je),
           label=as.integer(integer(nv)),
           status=as.integer(integer(1L)),
           PACKAGE="spatstat.geom")
  if(zz$status != 0)
    stop(paste("Internal error:", algoname, "did not converge"), call.=FALSE)
  return(zz$label)
}
  

# .................................................

is.connected <- function(X, ...) {
  UseMethod("is.connected")
}

is.connected.default <- function(X, ...) {
  y <- connected(X, ...)
  npieces <- length(levels(y))
  if(npieces == 0)
    stop("Unable to determine connectedness")
  return(npieces == 1)
}

is.connected.ppp <- function(X, R, ...) {
  lab <- connected(X, R, internal=TRUE)
  npieces <- length(unique(lab))
  return(npieces == 1)
}
