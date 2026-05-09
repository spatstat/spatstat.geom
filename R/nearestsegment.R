#'
#'  nearestsegment.R
#'
#'  $Revision: 1.15 $  $Date: 2026/05/07 05:50:13 $
#'
#' Given a point pattern X and a line segment pattern Y,
#' for each point x of X, determine which segment of Y is closest to x
#' and find the point on Y closest to x.
#'

nearestsegment <- function(X,Y) {
  return(ppllengine(X,Y,"identify"))
}

project2segment <- function(X, Y) {
  return(ppllengine(X,Y,"project"))
}
  
ppllengine <- function(X, Y, action="project",
                       method=c("C", "interpreted"),
                       check=FALSE) {
  stopifnot(is.ppp(X))
  stopifnot(is.psp(Y))
  action <- match.arg(action, c("distance", "identify", "project"))
  #' deal with empty patterns
  if(X$n == 0) {
    nowt <- numeric(0)
    none <- integer(0)
    switch(action,
           identify = return(none),
           distance = return(list(dist=nowt, which=none)),
           project  = return(list(Xproj=X, mapXY=none, d=nowt, tp=nowt)))
  }
  if(Y$n == 0)
    stop("Segment pattern Y contains 0 segments; projection undefined")
  #'              
  #' determine which segment lies closest to each point
  huge <- max(diameter(as.rectangle(as.owin(X))),
              diameter(as.rectangle(as.owin(Y))))
  XX <- as.matrix(as.data.frame(unmark(X)))
  YY <- as.matrix(as.data.frame(unmark(Y)))
  #' 
  switch(action,
         identify = {
           d <- distppllmin(XX, YY, huge^2)
           result <- d$min.which
         },
         distance = {
           d <- distppllmin(XX, YY, huge^2)
           result <- data.frame(dist=d$min.d, which=d$min.which)
         },
         project = {
           method <- match.arg(method)
           switch(method,
                  C = {
                    a <- NNdist2segments(XX[,1], XX[,2],
                                         YY[,1], YY[,2], YY[,3], YY[,4],
                                         bigvalue=huge^2,
                                         wantindex=TRUE, wantproj=TRUE)
                    Xproj <- ppp(a$xproj, a$yproj,
                                 window=X$window, marks=X$marks,
                                 check=check)
                    result <- list(Xproj=Xproj,
                                   mapXY=a$index,
                                   d=sqrt(a$dist2),
                                   tp=a$tproj)
                  },
                  interpreted = {
                    d <- distppllmin(XX, YY, huge^2)
                    mapXY <- d$min.which
                    #' combine relevant rows of data
                    alldata <- as.data.frame(cbind(XX, YY[mapXY, ,drop=FALSE]))
                    colnames(alldata) <- c("x", "y", "x0", "y0", "x1", "y1")
                    #' coordinate geometry
                    dx <- with(alldata, x1-x0)
                    dy <- with(alldata, y1-y0)
                    leng <- sqrt(dx^2 + dy^2)
                    #' rotation sines & cosines (may include 0/0)
                    co <- dx/leng
                    si <- dy/leng
                    #' vector to point from first endpoint of segment
                    xv <- with(alldata, x - x0)
                    yv <- with(alldata, y - y0)
                    #' rotate coordinate system so x axis is
                    #' parallel to line segment
                    xpr <- xv * co + yv * si
                    #' determine whether projection is
                    #' an endpoint or interior point
                    ok <- is.finite(xpr)
                    left <- !ok | (xpr <= 0)
                    right <- ok &  (xpr >= leng)
                    #' location of projected point in rotated coordinates
                    xr <- with(alldata,
                               ifelseAX(left, 0, ifelseXY(right, leng, xpr)))
                    #' back to standard coordinates
                    xproj <- with(alldata, x0 + ifelseXB(ok, xr * co, 0))
                    yproj <- with(alldata, y0 + ifelseXB(ok, xr * si, 0))
                    Xproj <- ppp(xproj, yproj, window=X$window, marks=X$marks,
                                 check=check)
                    #' parametric coordinates
                    tp <- xr/leng
                    tp[!is.finite(tp)] <- 0
                    result <- list(Xproj=Xproj, mapXY=mapXY, d=d$min.d, tp=tp)
                  })
           })
  return(result)
}

