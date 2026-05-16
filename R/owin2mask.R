#'
#'     owin2mask.R
#'
#'    Mask approximations which are guaranteed to be entirely inside
#'    or entirely covering the original window.
#'
#'    $Revision: 1.12 $  $Date: 2026/05/16 10:45:08 $
#'

owin2mask <- function(w,
                      ...,
                      op=c("sample", "notsample",
                           "cover", "inside", "uncover", "outside",
                           "boundary", "majority", "minority"),
                      W=w
                      ) {
  op <- if(!is.null(op)) match.arg(op) else "sample"
  if(!missing(w) && !missing(W))
    warning("Both arguments w and W were given in owin2mask", call.=FALSE)
  if(is.mask(W) && (length(list(...)) == 0)) {
    ## W is already a mask and there is no change to the raster
    switch(op,
           sample = ,
           cover = ,
           majority = ,
           inside = { return(W) },
           notsample = ,
           uncover = ,
           minority = ,
           outside = { return(complement.owin(W)) },
           boundary = { M <- W })
  } else {
    ## convert to mask, ignoring unrecognised arguments
    M <- AsMaskInternal(w=W, ...)
  }
  ## (M consists of all pixels whose centres are inside W)
  
  ## Do more processing
  switch(op,
         sample = ,
         notsample = {
           ## nothing
         },
         inside   = ,
         outside  = ,
         cover    = ,
         uncover  = ,
         boundary = {
           ## convert the boundary to a mask
           P <- as.polygonal(W)
           B <- psp2mask(edges(P), xy=M)
         },
         majority = ,
         minority = {
           ## compute the fraction of occupied area in each pixel
           U <- pixellate(W, M, DivideByPixelArea=TRUE)
         })

  ## Finally determine the mask
  R <- switch(op,
              sample    = M,
              notsample = complement.owin(M),
              inside   = setminus.owin(M, B, rescue=FALSE),
              outside  = setminus.owin(complement.owin(M), B, rescue=FALSE),
              cover    = union.owin(M, B, rescue=FALSE),
              uncover  = union.owin(complement.owin(M), B, rescue=FALSE),
              boundary = B,
              majority = levelset(U, 0.5, ">="), 
              minority = levelset(U, 0.5, "<"))
  return(R)
}

