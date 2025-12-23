#
#  is.subset.owin.R
#
#  $Revision: 1.18 $   $Date: 2025/12/23 01:45:47 $
#
#  Determine whether a window is a subset of another window
#
#  is.subset.owin()
#

is.subset.owin <- local({
  
  is.subset.owin <- function(A, B) {
    A <- as.owin(A)
    B <- as.owin(B)

    if(identical(A, B))
      return(TRUE)

    if(is.empty(A)) return(TRUE)
    if(is.empty(B)) return(FALSE)

    A <- rescue.rectangle(A)
    B <- rescue.rectangle(B)
  
    if(is.rectangle(B)) {
      # Some cases can be resolved using convexity of B
    
      # (1) A is also a rectangle
      if(is.rectangle(A)) {
        xx <- A$xrange[c(1L,2L,2L,1L)]
        yy <- A$yrange[c(1L,1L,2L,2L)]
        ok <- inside.owin(xx, yy, B)
        return(all(ok))
      } 
      # (2) A is polygonal
      # Then A is a subset of B iff,
      # for every constituent polygon of A with positive sign,
      # the vertices are all in B
      if(is.polygonal(A)) {
        ok <- unlist(lapply(A$bdry, okpolygon, B=B))
        return(all(ok))
      }
      # (3) Feeling lucky
      # Test whether the bounding box of A is a subset of B
      # Then a fortiori, A is a subset of B
      AA <- boundingbox(A)
      if(is.subset.owin(AA, B))
        return(TRUE)
    }

    if(!is.mask(A) && !is.mask(B)) {
      ## rectangles or polygonal domains
      if(!all(inside.owin(vertices(A), , B)))
        return(FALSE)
      ## all vertices of A are inside B.
      if(is.convex(B)) return(TRUE)
      ## check whether the boundaries are disjoint
      if(!anycrossing.psp(edges(A), edges(B))) {
        ## Disjoint boundary crossings sufficient if B has no holes
        if(length(B$bdry) == 1 || !any(sapply(B$bdry, is.hole.xypolygon)))
          return(TRUE)
        ## Compare area of intersection with area of putative subset
        ## (use '>=' instead of '==' because of numerical rounding error)
        areaA <- area(A)
        if(overlap.owin(A,B) >= areaA ||
           overlap.owin(B,A) >= areaA) return(TRUE)
      }
      ## continue...
    }
  
   # Discretise
    a <- as.mask(A)
    b <- as.mask(B)
    rxy <- rasterxy.mask(a, drop=TRUE)
    xx <- rxy$x
    yy <- rxy$y
    ok <- inside.owin(xx, yy, b)
    return(all(ok))
    
  }

  okpolygon <- function(a, B) {
    if(Area.xypolygon(a) < 0) return(TRUE)
    ok <- inside.owin(a$x, a$y, B)
    return(all(ok))
  }

  is.subset.owin
})
