#'
#'     owin2mask.R
#'
#'    Mask approximations which are guaranteed to be entirely inside
#'    or entirely covering the original window.
#'
#'    $Revision: 1.7 $  $Date: 2023/04/25 08:32:13 $
#'

owin2mask <- function(W, 
                      op=c("sample", "notsample",
                           "cover", "inside", "uncover", "outside",
                           "boundary", "majority", "minority"),
                      ...) {
  op <- match.arg(op)
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
    ## convert to mask
    M <- as.mask(W, ...)
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
           B <- as.mask.psp(edges(P), xy=M)
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
              inside   = setminus.owin(M, B),
              outside  = setminus.owin(complement.owin(M), B),
              cover    = union.owin(M, B),
              uncover  = union.owin(complement.owin(M), B),
              boundary = B,
              majority = levelset(U, 0.5, ">="), 
              minority = levelset(U, 0.5, "<"))
  return(R)
}

