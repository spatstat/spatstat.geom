#'
#'     owin2mask.R
#'
#'    Mask approximations which are guaranteed to be entirely inside
#'    or entirely covering the original window.
#'
#'    $Revision: 1.5 $  $Date: 2021/06/24 02:48:19 $
#'

owin2mask <- function(W, 
                      op=c("sample", "notsample",
                           "cover", "inside", "uncover", "outside",
                           "boundary"),
                      ...) {
  op <- match.arg(op)
  if(is.mask(W) && (length(list(...)) == 0)) {
    ## W is already a mask and there is no change to the raster
    switch(op,
           sample = ,
           cover = ,
           inside = { return(W) },
           notsample = ,
           uncover = ,
           outside = { return(complement.owin(W)) },
           boundary = { M <- W })
  } else {
    ## convert to mask
    M <- as.mask(W, ...)
  }
  ## M consists of all pixels whose centres are inside W
  if(op == "sample") return(M)
  if(op == "notsample") return(complement.owin(M))
  ## probe the boundary
  P <- as.polygonal(W)
  B <- as.mask.psp(edges(P), xy=M)
  R <- switch(op,
              inside   = setminus.owin(M, B),
              outside  = setminus.owin(complement.owin(M), B),
              cover    = union.owin(M, B),
              uncover  = union.owin(complement.owin(M), B),
              boundary = B,
              stop("Unrecognised option"))
  return(R)
}

