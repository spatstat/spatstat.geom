\name{spatstat.geom-internal} 
\title{Internal spatstat.geom functions}
\alias{[.diagramobj}
\alias{[.pp3}
\alias{[.splitppx}
\alias{[<-.splitppx}
\alias{acedist.show}
\alias{acedist.noshow}
\alias{affinexy}
\alias{affinexypolygon}
\alias{allElementsIdentical}
\alias{anycrossing.psp}
\alias{applytolayers}
\alias{applyPolyclipArgs}
\alias{areaGain.diri}       
\alias{areaGain.grid}       
\alias{areaLoss.diri}       
\alias{areaLoss.grid}
\alias{areaLoss.poly}
\alias{AsymmDistance.psp}
\alias{as.breakpts}
\alias{as.character.unitname}
\alias{as.data.frame.ppplist}
\alias{as.double.im}
\alias{as.imlist}
\alias{as.linimlist}
\alias{as.list.hyperframe}
\alias{as.listof}
\alias{as.ppplist}
\alias{as.unitname}
\alias{avenndist}
\alias{bbEngine}
\alias{bdry.mask}
\alias{boundingbox.list}
\alias{bounding.box3}
\alias{break.holes}
\alias{breakpts}
\alias{breakpts.from.r}
\alias{cartesian}
\alias{cellmiddles}
\alias{checkbigmatrix}          
\alias{checkfields}          
\alias{check.arc}
\alias{check.finespacing}
\alias{check.hist.lengths}
\alias{circunion}
\alias{clip.psp}
\alias{cliprect.psp}
\alias{clippoly.psp}
\alias{closethresh}
\alias{cocoEngine}
\alias{coerce.marks.numeric}
\alias{crosspairquad}
\alias{cobble.xy}
\alias{codetime}
\alias{col.args.to.grey}
\alias{commonPolyclipArgs}
\alias{conform.imagelist}
\alias{countingweights}
\alias{default.n.tiling}
\alias{default.ntile}
\alias{dimnames.hyperframe}
\alias{dimnames<-.hyperframe}
\alias{dflt.redraw}
\alias{diagramobj}
\alias{dim.hyperframe}               
\alias{dim.im}               
\alias{dim.owin}
\alias{dist2dpath}
\alias{do.as.im}
\alias{do.call.plotfun}
\alias{drawSignedPoly}
\alias{emptywindow}
\alias{equalpairs.quad}          
\alias{equals.quad}          
\alias{equalsfun.quad}          
\alias{erodemask}
\alias{even.breaks.owin}
\alias{exactdt}              
\alias{exactPdt}
\alias{existsSpatstatVariable}
\alias{expandSpecialLists}
\alias{fakemaintitle}
\alias{fft2D}
\alias{fftwAvailable}
\alias{fillNA}
\alias{flipxypolygon}
\alias{format.numberwithunit}
\alias{framebottomleft}
\alias{gammabreaks}
\alias{genericNNdistBy}
\alias{getfields}
\alias{getlastshift}
\alias{getSpatstatVariable}
\alias{gridindex}            
\alias{grid1index}
\alias{grow.mask}
\alias{handle.r.b.args}
\alias{head.hyperframe}
\alias{hsvNA}
\alias{IdenticalRowPair}
\alias{IdenticalRows}
\alias{idorempty}
\alias{imageOp}
\alias{imagelistOp}
\alias{inpoint}
\alias{inside.arc}
\alias{interpretAsOrigin}
\alias{intX.owin}     
\alias{intX.xypolygon}     
\alias{intY.owin}     
\alias{intY.xypolygon}
\alias{invokeColourmapRule}
\alias{is.col.argname}
\alias{is.data}
\alias{is.fv}
\alias{is.hyperframe}
\alias{is.imlist}
\alias{is.infline}
\alias{is.marked.default}    
\alias{is.marked.psp}    
\alias{is.marked.quad}
\alias{is.multitype.default}    
\alias{is.multitype.quad}    
\alias{is.ppplist}
\alias{is.pp3}
\alias{is.ppx}
\alias{is.psp}
\alias{is.quad}
\alias{is.sob}
\alias{is.solist}
\alias{is.tess}
\alias{is.vanilla}
\alias{levels.im}
\alias{levels<-.im}
\alias{levelsAsFactor}
\alias{listof}
\alias{logi.dummy}
\alias{markvaluetype}
\alias{packupNNdata}
\alias{parbreak}
\alias{perspVisible}
\alias{plan.legend.layout}
\alias{pointweights}
\alias{polytileareaEngine}
\alias{putSpatstatVariable}
\alias{lookup.im}
\alias{majorminorversion}
\alias{make.even.breaks}
\alias{makeunitname}
\alias{markappend}
\alias{markcbind}
\alias{markformat}
\alias{markformat.ppp}
\alias{markformat.ppx}
\alias{markformat.psp}
\alias{markformat.default}
\alias{mark.scale.default}
\alias{markspace.integral}
\alias{marks.default}           
\alias{marks.quad}           
\alias{\%mapp\%} %DoNotExport
%NAMESPACE export("%mapp%")
\alias{markappendop}
\alias{marksubset}
\alias{markreplicateop}
\alias{\%mrep\%} %DoNotExport
%NAMESPACE export("%mrep%")
\alias{marksubsetop}
\alias{\%msub\%} %DoNotExport
%NAMESPACE export("%msub%")
\alias{mask2df}
\alias{maxflow}
\alias{meanX.owin}
\alias{meanY.owin}
\alias{MinimalTess}
\alias{multiplicityNumeric}
\alias{multiply.only.finite.entries}
\alias{na.handle.im}
\alias{names.hyperframe}
\alias{names<-.hyperframe}
\alias{nearest.pixel}
\alias{nearest.valid.pixel}
\alias{n.quad}
\alias{numberwithunit}
\alias{numeric.columns}
\alias{onecolumn}
\alias{owinpolycheck}
\alias{owinpoly2mask}
\alias{owin2polypath}
\alias{param.quad}
\alias{PDtoNN}
\alias{pickoption}
\alias{plotEachLayer}
\alias{plot3Dpoints}
\alias{plotPolygonBdry}
\alias{plotWidthMap}
\alias{plot.barplotdata}
\alias{plot.indicfun}
\alias{ppllengine}
\alias{pppdist.mat}
\alias{pppdist.prohorov}
\alias{ppsubset}
\alias{prepareTitle}
\alias{print.anylist}    
\alias{print.colourmap}
\alias{print.distfun}
\alias{print.ewcdf}
\alias{print.funxy}       
\alias{print.hyperframe}
\alias{print.indicfun}       
\alias{print.layered}
\alias{print.lut}
\alias{print.metric}
\alias{print.metricfun}
\alias{print.nnfun}
\alias{print.numberwithunit}
\alias{print.onearrow}
\alias{print.pppmatching}
\alias{print.simplepanel}
\alias{print.solist}
\alias{print.splitppp}
\alias{print.splitppx}
\alias{print.summary.distfun}
\alias{print.summary.funxy}
\alias{print.summary.hyperframe}
\alias{print.summary.listof}
\alias{print.summary.logiquad}
\alias{print.summary.lut}
\alias{print.summary.owin}
\alias{print.summary.ppp}
\alias{print.summary.psp}
\alias{print.summary.solist}
\alias{print.summary.splitppp}
\alias{print.summary.splitppx}
\alias{print.summary.unitname}
\alias{print.symbolmap}       
\alias{print.textstring}
\alias{print.texturemap}
\alias{print.tess}
\alias{print.timed}
\alias{print.yardstick}
\alias{project3Dhom}
\alias{putlastshift}
\alias{qtPrepareCoordinate}
\alias{quad}
\alias{quantilefun.interpolatedCDF}
\alias{rasterfilter}
\alias{recognise.spatstat.type}
\alias{rectquadrat.breaks}
\alias{rectquadrat.countEngine}
\alias{remove.identical.pairs}
\alias{repair.image.xycoords}
\alias{resolve.stringsAsFactors}
\alias{rgbNA}
\alias{ruletextline}
\alias{quadscheme.replicated}
\alias{quadscheme.spatial}
\alias{pointgrid}
\alias{rastersample}
\alias{rasterx.mask}
\alias{rastery.mask}
\alias{rasterxy.mask}
\alias{rasterx.im}
\alias{rastery.im}
\alias{rasterxy.im}
\alias{rebound}
\alias{rebound.im}
\alias{rebound.ppp}
\alias{rebound.psp}
\alias{rebound.owin}
\alias{repair.old.factor.image}
\alias{restrict.mask}
\alias{reversePolyclipArgs}
\alias{rotxy}
\alias{rotxypolygon}
\alias{row.names.hyperframe}
\alias{row.names<-.hyperframe}
\alias{safedeldir}
\alias{safeDevCapabilities}
\alias{safelookup}
\alias{scalardilate.breakpts}
\alias{scalardilate.diagramobj}
\alias{shift.diagramobj}              
\alias{shift.quadratcount}              
\alias{shiftxy}              
\alias{shiftxypolygon}              
\alias{simulationresult}
\alias{smudge}
\alias{sort.im}
\alias{spatstatDiagnostic}
\alias{spatstat.deldir.setopt}
\alias{spatstat.xy.coords}
\alias{store.versionstring.spatstat}
\alias{str.hyperframe}
\alias{summary.hyperframe}
\alias{summary.logiquad}
\alias{summary.lut}
\alias{summary.metric}
\alias{summary.pppmatching}
\alias{summary.ppx}
\alias{summary.splitppx}
\alias{superimposeMarks}
\alias{symbolmapdomain}
\alias{symbolmapparnames}
\alias{symbolmaptype}
\alias{tail.hyperframe}
\alias{thickSegments}
\alias{tilecentroids}        
\alias{trianglediameters}
\alias{trim.mask}
\alias{tweak.closepairs}
\alias{\%unit\%} %DoNotExport
%NAMESPACE export("%unit%")
\alias{unitname.default}
\alias{unitname<-.default}
\alias{unstackFilter}
\alias{update.im}
\alias{validradius}
\alias{validate.mask}        
\alias{validate.quad}        
\alias{versioncurrency.spatstat}
\alias{versionstring.spatstat}
\alias{verifyclass}
\alias{veryunique}
\alias{warn.once}
\alias{warn.no.metric.support}
\alias{warn.unsupported.args}
\alias{waxlyrical}
\alias{wrangle2image}
\alias{w.quad}               
\alias{x.quad}
\alias{y.quad}
\alias{xy.grid}
\alias{xtfrm.im}
\alias{XDtoNN}
\alias{xypolygon2psp}
\alias{xypolyselfint}
\alias{ZeroValue}
\alias{ZeroValue.im}




%%%%%%%
\description{
  Internal spatstat.geom functions.
}
\usage{
\method{[}{diagramobj}(x, \dots)
\method{[}{pp3}(x, i, drop, \dots)
\method{[}{splitppx}(x, \dots)
\method{[}{splitppx}(x, \dots) <- value
acedist.show(X, Y, n, d, timelag)
acedist.noshow(X, Y, n, d)
affinexy(X, mat, vec, invert)
affinexypolygon(p, mat, vec, detmat)
allElementsIdentical(x, entry)
anycrossing.psp(A,B)
applytolayers(L, FUN, \dots)
applyPolyclipArgs(x, p)
areaGain.diri(u, X, r, \dots, W, verbose)
areaGain.grid(u, X, r, \dots, W, ngrid)
areaLoss.diri(X, r, \dots, W, subset)
areaLoss.grid(X, r, \dots, W, subset,
                         method = c("count", "distmap"),
                         ngrid = spatstat.options("ngrid.disc"),
                         exact = FALSE)
areaLoss.poly(X, r, \dots, W, subset, splitem)
AsymmDistance.psp(X, Y, metric, method)
as.breakpts(\dots)
\method{as.character}{unitname}(x, \dots)
\method{as.data.frame}{ppplist}(x, row.names, \dots)
\method{as.double}{im}(x, \dots)
as.imlist(x, check)
as.linimlist(x, check)
\method{as.list}{hyperframe}(x, \dots)
as.listof(x)
as.ppplist(x, check)
as.unitname(s)
avenndist(X)
bbEngine(\dots)
bdry.mask(W)
\method{boundingbox}{list}(\dots)
bounding.box3(\dots)
break.holes(x, splitby, depth, maxdepth)
breakpts(val, maxi, even = FALSE, npos = NULL, step = NULL)
breakpts.from.r(r)
cartesian(pp, markset, fac = TRUE)
cellmiddles(W, nx, ny, npix, distances)
checkbigmatrix(n, m, fatal, silent)
checkfields(X,L)
check.arc(arc, fatal)
check.finespacing(r, eps, win, rmaxdefault, context, action, rname)
check.hist.lengths(hist,breaks)
circunion(arcs)
clip.psp(x, window, check, fragments)
cliprect.psp(x, window, fragments)
clippoly.psp(s, window, fragments)
closethresh(X,R,S,twice,\dots)
cocoEngine(nv, ie, je, algoname)
coerce.marks.numeric(X, warn)
crosspairquad(Q,rmax,what)
cobble.xy(x, y, f, fatal, \dots)
codetime(x, hms, what)
col.args.to.grey(x, \dots)
commonPolyclipArgs(\dots, p)
conform.imagelist(X, Zlist)
countingweights(id, areas, check = TRUE)
default.n.tiling(X, nd, ntile, npix, eps, random, quasi, verbose)
default.ntile(X)
\method{dimnames}{hyperframe}(x)
\method{dimnames}{hyperframe}(x) <- value
dflt.redraw(button, name, env)
diagramobj(X, \dots)
\method{dim}{hyperframe}(x)
\method{dim}{im}(x)
\method{dim}{owin}(x)
dist2dpath(dist, method="C")
do.as.im(x, action, \dots, W, eps, dimyx, xy, rule.eps, na.replace)
do.call.plotfun(fun, arglist, \dots, envir)
drawSignedPoly(x,y,pars,sgn)
emptywindow(w)
equalpairs.quad(Q)
equals.quad(Q)          
equalsfun.quad(Q)
erodemask(w,r,strict)
even.breaks.owin(w)
exactdt(X, \dots)              
exactPdt(w)
existsSpatstatVariable(name)
expandSpecialLists(x, special)
fakemaintitle(bb, main, \dots)
fft2D(z, inverse, west)
fftwAvailable()
fillNA(x, value)
flipxypolygon(p)
\method{format}{numberwithunit}(x, \dots, collapse, modifier)
framebottomleft(w)
gammabreaks(ra, n, gamma)
genericNNdistBy(X, by, k)
getfields(X, L, fatal = TRUE)
getlastshift(X)
getSpatstatVariable(name, default)
gridindex(x, y, xrange, yrange, nx, ny)            
grid1index(x, xrange, nx)
grow.mask(M, xmargin=0, ymargin=xmargin)
handle.r.b.args(r = NULL, breaks = NULL, window, pixeps = NULL, rmaxdefault)
\method{head}{hyperframe}(x,n,\dots)
hsvNA(h, s, v, alpha)
IdenticalRowPair(i,j,a,b)
IdenticalRows(i,j,a,b)
idorempty(w, r, caller)
imageOp(e1, e2, op)
imagelistOp(e1, e2, op)
inpoint(W)
inside.arc(theta, arc)
interpretAsOrigin(x, W)
intX.owin(w)
intX.xypolygon(polly)
intY.owin(w)
intY.xypolygon(polly)
invokeColourmapRule(colfun, x, \dots, zlim, colargs)
is.col.argname(x)
is.data(Q)
is.fv(x)
is.hyperframe(x)
is.imlist(x)
is.infline(x)
\method{is.marked}{default}(\dots)  
\method{is.marked}{psp}(X, \dots)
\method{is.marked}{quad}(X, na.action="warn", \dots)
\method{is.multitype}{default}(X, \dots)  
\method{is.multitype}{quad}(X, na.action="warn", \dots)
is.ppplist(x)
is.pp3(x)
is.ppx(x)
is.psp(x)
is.quad(x)
is.solist(x)
is.sob(x)
is.tess(x)
is.vanilla(u)
\method{levels}{im}(x)
\method{levels}{im}(x) <- value
levelsAsFactor(x)
listof(\dots)
logi.dummy(X, dummytype, nd, mark.repeat, \dots)
markvaluetype(x)
packupNNdata(NND, NNW, what, k)
parbreak(terse)
perspVisible(x, y, z, M)
plan.legend.layout(B, \dots, side, sep, size, sep.frac, size.frac,
                   started, map)
pointweights(X, \dots, weights, parent)
polytileareaEngine(P, xrange, yrange, nx, ny, DivideByPixelArea)
putSpatstatVariable(name, value)
lookup.im(Z, x, y, naok, strict)
majorminorversion(v)
make.even.breaks(bmax, npos, bstep)
makeunitname(sing, plur, mul)
markappend(\dots)
markcbind(\dots)
markformat(x)
\method{markformat}{ppp}(x) 
\method{markformat}{ppx}(x) 
\method{markformat}{psp}(x) 
\method{markformat}{default}(x) 
mark.scale.default(marx, w, \dots, markscale, maxsize, meansize, characters)
markspace.integral(X)
\method{marks}{default}(x, \dots)
\method{marks}{quad}(x, dfok=FALSE, \dots)
markappendop(x, y)
x \%mapp\% y
marksubset(x, index, format)
marksubsetop(x, i)
x \%msub\% i
markreplicateop(x, n)
x \%mrep\% n
mask2df(w)
maxflow(costm)
meanX.owin(w)            
meanY.owin(w)
MinimalTess(W, \dots)
multiplicityNumeric(x)
multiply.only.finite.entries(x, a)
na.handle.im(X, na.replace)
\method{names}{hyperframe}(x)
\method{names}{hyperframe}(x) <- value
nearest.pixel(x, y, Z)
nearest.valid.pixel(x, y, Z, method, nsearch)
n.quad(Q)
numberwithunit(x, u)
numeric.columns(M, logical, others)
onecolumn(m)
owinpolycheck(W, verbose=TRUE)
owinpoly2mask(w, rasta, check=TRUE)
owin2polypath(w)
param.quad(Q)
PDtoNN(d, what, k, \dots)
pickoption(what="option", key, keymap, \dots,
           exact=FALSE, list.on.err=TRUE, die=TRUE, multi=FALSE,
           allow.all=TRUE)
plotEachLayer(x, \dots, main, plotargs, add, show.all, do.plot)
plot3Dpoints(xyz, eye, org,
             \dots,
             type, xlim, ylim, zlim,
             add, box, main, cex, box.back, box.front)
plotPolygonBdry(x, \dots)
plotWidthMap(bb.leg, zlim, phys.scale, leg.scale, leg.side, leg.args, grafpar)
\method{plot}{barplotdata}(x, \dots)
\method{plot}{indicfun}(x, W, \dots, main)
ppllengine(X, Y, action="project", check=FALSE)
pppdist.mat(X, Y, cutoff = 1, q = 1, matching = TRUE,
            precision = 9, approximation = 10)
pppdist.prohorov(X, Y, n, dfix, type, cutoff, matching,
            ccode, auction, precision, approximation) 
ppsubset(X, I, Iname, fatal)
prepareTitle(main)
\method{print}{anylist}(x, \dots)
\method{print}{colourmap}(x, \dots)
\method{print}{distfun}(x, \dots)
\method{print}{ewcdf}(x, digits, \dots)
\method{print}{funxy}(x, \dots)
\method{print}{hyperframe}(x, \dots)
\method{print}{indicfun}(x, \dots)
\method{print}{layered}(x, \dots)
\method{print}{lut}(x, \dots)
\method{print}{metric}(x, \dots)
\method{print}{metricfun}(x, \dots)
\method{print}{nnfun}(x, \dots)
\method{print}{numberwithunit}(x, \dots)
\method{print}{onearrow}(x, \dots)
\method{print}{pppmatching}(x, \dots)
\method{print}{simplepanel}(x, \dots)
\method{print}{solist}(x, \dots)
\method{print}{splitppp}(x, \dots)
\method{print}{splitppx}(x, \dots)
\method{print}{summary.distfun}(x, \dots)
\method{print}{summary.funxy}(x, \dots)
\method{print}{summary.hyperframe}(x, \dots)
\method{print}{summary.listof}(x, \dots)
\method{print}{summary.logiquad}(x, \dots, dp=3)
\method{print}{summary.lut}(x, \dots)
\method{print}{summary.owin}(x, \dots)
\method{print}{summary.ppp}(x, \dots, dp)
\method{print}{summary.psp}(x, \dots)
\method{print}{summary.splitppp}(x, \dots)
\method{print}{summary.solist}(x, \dots)
\method{print}{summary.splitppx}(x, \dots)
\method{print}{summary.unitname}(x, \dots)
\method{print}{symbolmap}(x, \dots)
\method{print}{textstring}(x, \dots)
\method{print}{texturemap}(x, \dots)
\method{print}{tess}(x, \dots, brief=FALSE)
\method{print}{timed}(x, \dots)
\method{print}{yardstick}(x, \dots)
project3Dhom(xyz, eye, org, vert)
putlastshift(X, vec)
qtPrepareCoordinate(covname, W, origin)
quad(data, dummy, w, param)
\method{quantilefun}{interpolatedCDF}(x, \dots, type)
rasterfilter(X, f)
recognise.spatstat.type(x)
rectquadrat.breaks(xr, yr, nx = 5, ny = nx, xbreaks = NULL, ybreaks = NULL)
rectquadrat.countEngine(x, y, xbreaks, ybreaks, weights)
remove.identical.pairs(cl, imap, jmap)
repair.image.xycoords(x)
resolve.stringsAsFactors(stringsAsFactors)
rgbNA(red, green, blue, alpha, maxColorValue)
ruletextline(ch, n, terse)
quadscheme.replicated(data, dummy, method, \dots)
quadscheme.spatial(data, dummy, method, \dots)
pointgrid(W, ngrid)
rastersample(X, Y)
rasterx.mask(w, drop)
rastery.mask(w, drop)
rasterxy.mask(w, drop)
rasterx.im(x)
rastery.im(x)
rasterxy.im(x, drop)
rebound(x, rect)
\method{rebound}{im}(x, rect)  
\method{rebound}{ppp}(x, rect) 
\method{rebound}{psp}(x, rect) 
\method{rebound}{owin}(x, rect)
repair.old.factor.image(x)
restrict.mask(M, W)
reversePolyclipArgs(x, p)
rotxy(X, angle = pi/2)
rotxypolygon(p, angle = pi/2)
\method{row.names}{hyperframe}(x)
\method{row.names}{hyperframe}(x) <- value
safedeldir(X)
safeDevCapabilities()
safelookup(Z, x, factor, warn)
\method{scalardilate}{breakpts}(X, f, \dots)
\method{scalardilate}{diagramobj}(X, f, \dots)
\method{shift}{diagramobj}(X, \dots)
\method{shift}{quadratcount}(X, \dots)
shiftxy(X, vec = c(0, 0))
shiftxypolygon(p, vec = c(0, 0))
simulationresult(resultlist, nsim, drop, NameBase)
smudge(X)
\method{sort}{im}(x, \dots)
spatstatDiagnostic(msg)
spatstat.deldir.setopt(use.trigrafS, use.trigraf, debug.delaunay)
spatstat.xy.coords(x, y)
store.versionstring.spatstat()
\method{str}{hyperframe}(object, \dots)
\method{summary}{hyperframe}(object, \dots, brief=FALSE)
\method{summary}{logiquad}(object, \dots, checkdup=FALSE)
\method{summary}{lut}(object, \dots)
\method{summary}{metric}(object, \dots)
\method{summary}{pppmatching}(object, \dots)
\method{summary}{ppx}(object, \dots)
\method{summary}{splitppx}(object, \dots)
superimposeMarks(arglist, nobj)
symbolmapdomain(x)
symbolmapparnames(x)
symbolmaptype(x)
\method{tail}{hyperframe}(x,n,\dots)
thickSegments(x, widths, \dots,
             add, main, do.plot, show.all, show.window,
             scale, adjust, negative.args,
             legend, leg.side, leg.sep, leg.wid, leg.args, leg.scale,
             zlim, box)
tilecentroids(W, nx, ny)
trianglediameters(iedge, jedge, edgelength, \dots,
                  nvert, dmax, check)
trim.mask(M, R, tolerant)
tweak.closepairs(cl, rmax, i, deltax, deltay, deltaz)
x \%unit\% u
\method{unitname}{default}(x) 
\method{unitname}{default}(x) <- value 
unstackFilter(x)
\method{update}{im}(object, \dots) 
validradius(r, caller)
validate.mask(w, fatal=TRUE)        
validate.quad(Q, fatal, repair, announce)
versioncurrency.spatstat(today, checkR)
versionstring.spatstat()
veryunique(z)
verifyclass(X, C, N = deparse(substitute(X)), fatal = TRUE)
warn.once(key, \dots)
warn.no.metric.support(caller, \dots, metric)
warn.unsupported.args(unsup, \dots)
waxlyrical(type, terse)
wrangle2image(values, template)
w.quad(Q)               
x.quad(Q)
y.quad(Q)
xy.grid(xr, yr, nx, ny, dx, dy)
\method{xtfrm}{im}(x)
xypolyselfint(p, eps, proper, yesorno, checkinternal)
XDtoNN(d, what, iX, iY, k, \dots)
xypolygon2psp(p, w, check)
ZeroValue(x)
\method{ZeroValue}{im}(x)





}
\details{
  These internal \pkg{spatstat.geom} functions should not be called
  directly by the user. Their names and capabilities may change
  without warning from one version of \pkg{spatstat.geom} to the next.
}
\value{
The return values of these functions are not documented,
and may change without warning.
}
\keyword{internal}
