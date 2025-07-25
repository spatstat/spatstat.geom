#
# clickpoly.R
#
#
# $Revision: 1.13 $  $Date: 2025/07/25 06:22:09 $
#
#

clickpoly <- function(add=FALSE, nv=NULL, np=1, ...,
                      snap.step=NULL, snap.origin=c(0,0)) {
  if((!add) | dev.cur() == 1L) {
    plot(0,0,type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), asp=1.0,
         axes=FALSE)
    rect(0,0,1,1)
  }
  spatstatLocator(0) ## check locator is enabled
  if(snappy <- !is.null(snap.step)) {
    ## draw guide lines
    usr <- par('usr')
    snap.step <- ensure2vector(snap.step)
    botleft <- snapxy(usr[1], usr[3], snap.step, snap.origin)
    upright <- snapxy(usr[2], usr[4], snap.step, snap.origin)
    abline(v=seq(botleft$x, upright$x, by=snap.step[1]), lty=3, col=8)
    abline(h=seq(botleft$y, upright$y, by=snap.step[2]), lty=3, col=8)
  }
  gon <- list()
  stopifnot(np >= 1)
  #
  for(i in 1:np) {
    if(np > 1)
      cat(paste(".... Polygon number", i, ".....\n"))
    if(!is.null(nv)) 
      cat(paste("click", nv, "times in window\n"))
    else
      cat(paste("to add points: click left mouse button in window\n",
                "      to exit: press ESC or click middle mouse button\n",
                "[The last point should NOT repeat the first point]\n"))
    xy <- do.call(spatstatLocator,
                  resolve.defaults(if(!is.null(nv)) list(n=nv) else list(),
                                   list(...),
                                   list(snap.step=snap.step,
                                        snap.origin=snap.origin),
                                   list(type="o")))
    if(Area.xypolygon(xy) < 0)
      xy <- lapply(xy, rev)
    gon[[i]] <- xy
    plotPolygonBdry(owin(poly=xy), ...)
  }
  result <- owin(poly=gon)
  plotPolygonBdry(result, ...)
  return(result)
}

clickbox <- function(add=TRUE, ...) {
  spatstatLocator(0) # check locator enabled
  cat("Click two corners of a box\n")
  if(!add) plot(owinInternalRect(), main="Click two corners of a box") 
  a <- try(spatstatLocator(1), silent=TRUE)
  if(inherits(a, "try-error")) {
    ## add=TRUE but there is no current plot
    plot.new()
    a <- spatstatLocator(1, ...)
  }
  abline(v=a$x)
  abline(h=a$y)
  b <- spatstatLocator(1, ...)
  abline(v=b$x)
  abline(h=b$y)
  ab <- concatxy(a, b)
  result <- owinInternalRect(range(ab$x), range(ab$y))
  plotPolygonBdry(result, ...)
  return(result)
}

plotPolygonBdry <- function(x, ...) {
  # filter appropriate arguments
  argh <- list(...)
  polyPars <- union(graphicsPars("lines"), graphicsPars("owin"))
  polyargs <- argh[names(argh) %in% polyPars]
  # change 'col' to 'border'
  nama <- names(polyargs)
  if(any(nama == "col") && !any(nama == "border"))
    names(polyargs)[nama == "col"] <- "border"
  # plot
  do.call(plot.owin,
          append(list(x=x, add=TRUE), polyargs))
}
