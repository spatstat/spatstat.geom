#'
#'     owin2mask.R
#'
#'    Mask approximations which are guaranteed to be entirely inside
#'    or entirely covering the original window.
#'
#'    $Revision: 1.14 $  $Date: 2026/05/17 02:47:11 $
#'

owin2mask <- function(w,
                      ...,
                      rule.pix = c("sample", "notsample",
                                   "cover", "inside", "uncover", "outside",
                                   "boundary", "majority", "minority"),
                      W=w) {
  if(missing(w)) w <- NULL
  W <- w %orifnull% W
  rule.pix <- if(is.null(rule.pix)) "sample" else match.arg(rule.pix)
  if(is.mask(W) && (length(list(...)) == 0)) {
    ## W is already a mask and there is no change to the raster
    switch(rule.pix,
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
  switch(rule.pix,
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
  R <- switch(rule.pix,
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

