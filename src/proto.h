#include <R.h>
#include <Rinternals.h>

/*
  Prototype declarations for all native routines in spatstat.core package

  Automatically generated - do not edit! 

*/

/*
  
                  Functions invoked by .C

*/

void areadifs(double *, int *, double *, double *, int *, int *, double *); 
void areaBdif(double *, int *, double *, double *, int *, int *, double *, double *, double *, double *, double *);
void xysegint(int *, double *, double *, double *, double *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, int *);
void Fclosepairs(int *, double *, double *, double *, int *, int *, int *, int *, double *, double *, double *, double *, double *, double *, double *, int *); 
void paircount(int *, double *, double *, double *, int *); 
void Fclosepairs(int *, double *, double *, double *, int *, int *, int *, int *, double *, double *, double *, double *, double *, double *, double *, int *); 
void crosscount(int *, double *, double *, int *, double *, double *, double *, int *); 
void Fcrosspairs(int *, double *, double *, int *, double *, double *, double *, int *, int *, int *, int *, double *, double *, double *, double *, double *, double *, double *, int *);
void coco8dbl(double *, int *, int *); void coco4dbl(double *, int *, int *); 
void coco8int(int *, int *, int *); void coco4int(int *, int *, int *); 
void cocoGraph(int *, int *, int *, int *, int *, int *);
void mdtPconv(double *, double *, double *, double *, int *, int *, int *, int *, int *, int *, double *, double *, int *, double *, int *, int *); 
void mdtPconv(double *, double *, double *, double *, int *, int *, int *, int *, int *, int *, double *, double *, int *, double *, int *, int *);
void trigrafS(int *, int *, int *, int *, int *, int *, int *, int *, int *, int *); 
void trigraf(int *, int *, int *, int *, int *, int *, int *, int *, int *, int *); 
void Idist2dpath(int *, int *, int *, int *, int *, int *, int *);
void discareapoly(int *, double *, double *, int *, double *, int *, double *, double *, double *, double *, double *, double *);
void Ddist2dpath(int *, double *, int *, double *, double *, int *, int *);
void D3pairdist(int *, double *, double *, double *, int *, double *); 
void D3pairPdist(int *, double *, double *, double *, double *, double *, double *, int *, double *); 
void nnd3D(int *, double *, double *, double *, double *, int *, double *); 
void knnd3D(int *, int *, double *, double *, double *, double *, int *, double *); 
void nnw3D(int *, double *, double *, double *, double *, int *, double *); 
void knnw3D(int *, int *, double *, double *, double *, double *, int *, double *); 
void D3crossdist(int *, double *, double *, double *, int *, double *, double *, double *, int *, double *); 
void D3crossPdist(int *, double *, double *, double *, int *, double *, double *, double *, double *, double *, double *, int *, double *);
void Cpairdist(int *, double *, double *, int *, double *); 
void CpairPdist(int *, double *, double *, double *, double *, int *, double *); 
void Ccrossdist(int *, double *, double *, int *, double *, double *, int *, double *); 
void CcrossPdist(int *, double *, double *, int *, double *, double *, double *, double *, int *, double *);
void nndMD(int *, int *, double *, double *, double *); 
void knndMD(int *, int *, int *, double *, double *, double *); 
void nnwMD(int *, int *, double *, double *, int *, double *); 
void knnwMD(int *, int *, int *, double *, double *, int *, double *); 
void nnXwMD(int *, int *, double *, int *, double *, double *, int *, double *); 
void nnXxMD(int *, int *, double *, int *, int *, double *, int *, double *, int *, double *); 
void knnXwMD(int *, int *, double *, int *, double *, int *, double *, int *, double *); 
void knnXxMD(int *, int *, double *, int *, int *, double *, int *, int *, double *, int *, double *);
void distmapbin(int *, double *, double *, double *, double *, int *, int *, int *, double *, double *);
void tabsumweight(int *, double *, double *, int *, double *, double *);
void exact_dt_R(double *, double *, int *, double *, double *, double *, double *, int *, int *, int *, int *, double *, int *, double *);
void ps_exact_dt_R(double *, double *, double *, double *, int *, int *, int *, int *, int *, double *, int *, int *, double *);
void fardist2grid(int *, double *, double *, int *, double *, double *, int *, double *, double *, double *); 
void fardistgrid(int *, double *, double *, int *, double *, double *, int *, double *, double *, double *);
void hasXclose(int *, double *, double *, double *, int *); 
void hasXpclose(int *, double *, double *, double *, double *, int *); 
void hasXYclose(int *, double *, double *, int *, double *, double *, double *, int *); 
void hasXYpclose(int *, double *, double *, int *, double *, double *, double *, double *, int *); 
void hasX3close(int *, double *, double *, double *, double *, int *); 
void hasX3pclose(int *, double *, double *, double *, double *, double *, int *); 
void hasXY3close(int *, double *, double *, double *, int *, double *, double *, double *, double *, int *); 
void hasXY3pclose(int *, double *, double *, double *, int *, double *, double *, double *, double *, double *, int *);
void hotrodInsul(int *, double *, double *, double *, double *, int *, double *); 
void hotrodAbsorb(int *, double *, double *, double *, double *, int *, double *);
void nearestvalidpixel(int *, double *, double *, int *, int *, double *, int *, int *, int *, int *);
void mdtPOrect(double *, double *, double *, double *, int *, int *, int *, int *, int *, double *, int *, double *, int *, int *); 
void mdtPOrect(double *, double *, double *, double *, int *, int *, int *, int *, int *, double *, int *, double *, int *, int *);
void minPnnd2(int *, double *, double *, double *, double *); 
void minnnd2(int *, double *, double *, double *, double *); 
void maxPnnd2(int *, double *, double *, double *, double *); 
void maxnnd2(int *, double *, double *, double *, double *);
void nnX3Dinterface(int *, double *, double *, double *, int *, int *, double *, double *, double *, int *, int *, int *, int *, double *, int *, double *); 
void knnX3Dinterface(int *, double *, double *, double *, int *, int *, double *, double *, double *, int *, int *, int *, int *, int *, double *, int *, double *);
void nnXinterface(int *, double *, double *, int *, int *, double *, double *, int *, int *, int *, int *, double *, int *, double *); 
void knnXinterface(int *, double *, double *, int *, int *, double *, double *, int *, int *, int *, int *, int *, double *, int *, double *);
void nndistsort(int *, double *, double *, double *, double *); 
void knndsort(int *, int *, double *, double *, double *, double *); 
void nnwhichsort(int *, double *, double *, int *, double *); 
void knnwhich(int *, int *, double *, double *, int *, double *);
void nnGinterface(int *, double *, double *, int *, double *, double *, int *, double *, double *, int *, int *, double *, int *, double *); 
void knnGinterface(int *, double *, double *, int *, double *, double *, int *, double *, double *, int *, int *, int *, double *, int *, double *);
void poly2imA(int *, int *, double *, double *, int *, double *, int *);
void xypsi(int *, double *, double *, double *, double *, double *, double *, double *, int *, int *); 
void Cxypolyselfint(int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, int *);
void auctionbf(int *, int *, int *, double *, double *, int *, double *); 
void dwpure(int *, int *, int *, int *, int *, int *); 
void auctionbf(int *, int *, int *, double *, double *, int *, double *); 
void dwpure(int *, int *, int *, int *, int *, int *); 
void dinfty_R(int *, int *, int *); 
void dwpure(int *, int *, int *, int *, int *, int *); 
void dwpure(int *, int *, int *, int *, int *, int *);
void seg2pixI(int *, double *, double *, double *, double *, int *, int *, int *); 
void seg2pixL(int *, double *, double *, double *, double *, double *, double *, double *, int *, int *, double *); 
void seg2pixN(int *, double *, double *, double *, double *, double *, int *, int *, double *);
void xysegint(int *, double *, double *, double *, double *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, int *); 
void xysi(int *, double *, double *, double *, double *, int *, double *, double *, double *, double *, double *, int *); 
void xysiANY(int *, double *, double *, double *, double *, int *, double *, double *, double *, double *, double *, int *); 
void xysegXint(int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, int *); 
void xysxi(int *, double *, double *, double *, double *, double *, int *);
void Corput(int *, int *, double *);
void raster3filter(int *, int *, double *, double *, double *);
void uniqmapxy(int *, double *, double *, int *); 
void uniqmap2M(int *, double *, double *, int *, int *);
void anydupxy(int *, double *, double *, int *);
void poly2imI(double *, double *, int *, int *, int *, int *); 
void bdrymask(int *, int *, int *, int *); 
void discs2grid(int *, double *, double *, int *, double *, double *, int *, double *, double *, double *, int *);
/*

             Functions invoked by .Call

*/
SEXP close3pairs(SEXP, SEXP, SEXP, SEXP, SEXP); 
SEXP close3IJpairs(SEXP, SEXP, SEXP, SEXP, SEXP); 
SEXP close3IJDpairs(SEXP, SEXP, SEXP, SEXP, SEXP); 
SEXP cross3pairs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP); 
SEXP cross3IJpairs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP); 
SEXP cross3IJDpairs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP closePpair(SEXP, SEXP, SEXP, SEXP, SEXP); 
SEXP Vclosepairs(SEXP, SEXP, SEXP, SEXP); 
SEXP VcloseIJpairs(SEXP, SEXP, SEXP, SEXP); 
SEXP VcloseIJDpairs(SEXP, SEXP, SEXP, SEXP); 
SEXP altVclosepairs(SEXP, SEXP, SEXP, SEXP); 
SEXP altVcloseIJpairs(SEXP, SEXP, SEXP, SEXP); 
SEXP altVcloseIJDpairs(SEXP, SEXP, SEXP, SEXP); 
SEXP crossPpair(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP); 
SEXP Vcrosspairs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP); 
SEXP VcrossIJpairs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP); 
SEXP VcrossIJDpairs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP); 
SEXP Vclosethresh(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP trioxgraph(SEXP, SEXP, SEXP, SEXP); SEXP triograph(SEXP, SEXP, SEXP); 
SEXP trigraph(SEXP, SEXP, SEXP); SEXP triDgraph(SEXP, SEXP, SEXP, SEXP); 
SEXP triDRgraph(SEXP, SEXP, SEXP, SEXP, SEXP); 
SEXP graphVees(SEXP, SEXP, SEXP);
SEXP Cxysegint(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP); 
SEXP CxysegXint(SEXP, SEXP, SEXP, SEXP, SEXP); 
SEXP CxysegXint(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Cwhist(SEXP, SEXP, SEXP);
