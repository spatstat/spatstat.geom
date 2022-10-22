/*
  nearestpix.c

  Find the nearest TRUE pixel to a given (x,y) location

  Query locations (x,y) are transformed to coordinates
  in which the pixels are unit squares 
  and the pixel centres start at (0,0)

  $Revision: 1.3 $ $Date: 2022/10/21 10:43:01 $

 */

#include <R.h>
#include <R_ext/Utils.h>
#include <Rmath.h>
#include <math.h>

#include "chunkloop.h"
#include "yesno.h"

double sqrt(), fround();

void 
nearestvalidpixel(
  int *n,         /* number of query points */
  double *x,
  double *y,      /* coordinates of query points (transformed) */
  int *nr,
  int *nc,        /* matrix dimensions */
  double *aspect, /* aspect ratio (y/x) of original pixels */
  int *z,         /* entries of logical matrix */
  int *nsearch,   /* maximum permitted number of pixel steps on each axis */
  /* OUTPUTS */
  int *rr,
  int *cc         /* row and column indices (-1) of nearest pixel centre */
) {
  int maxchunk, N, Nrow, Ncol, maxrow, maxcol, maxsearch;
  double asp, xi, yi, ddd, ddi, huge, deltax, deltay;
  int i, row, col, zvalue;
  int rrr, ccc, rri, cci, startrow, endrow, startcol, endcol;
  
  N = *n;
  Nrow = *nr;
  Ncol = *nc;
  maxsearch = *nsearch;
  asp  = *aspect;
  
  maxrow = Nrow - 1;
  maxcol = Ncol - 1;
  huge = sqrt(((double) Ncol) * ((double) Ncol) +
	      asp * asp * ((double) Nrow) * ((double) Nrow));
  
  OUTERCHUNKLOOP(i, N, maxchunk, 8196) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, N, maxchunk, 8196) {
      xi = x[i];
      yi = y[i];
      row = (int) fround(yi, (double) 0);
      col = (int) fround(xi, (double) 0);
      if(row < 0) row = 0; else if(row > maxrow) row = maxrow;
      if(col < 0) col = 0; else if(col > maxcol) col = maxcol;
      zvalue = z[row + Nrow * col];
      if(zvalue != 0) {
	/* pixel is TRUE */
	rr[i] = row;
	cc[i] = col;
      } else {
	/* initialise result to NA */
	rri = cci = -1;
	ddi = huge;
	/* search neighbouring pixels */
	startrow = imax2(row - maxsearch, 0);
	endrow   = imin2(row + maxsearch, maxrow);
	startcol = imax2(col - maxsearch, 0);
	endcol   = imin2(col + maxsearch, maxcol);
	if(startrow <= endrow && startcol <= endcol) {
	  for(rrr = startrow; rrr <= endrow; rrr++) {
	    for(ccc = startcol; ccc <= endcol; ccc++) {
	      zvalue = z[rrr + Nrow * ccc];
	      if(zvalue != 0) {
		/* pixel is TRUE */
		deltax = xi - (double) ccc;
		deltay = asp * (yi - (double) rrr);
		ddd = sqrt(deltax * deltax + deltay * deltay);
		if(ddd < ddi) {
		  /* pixel is closer */
		  rri = rrr;
		  cci = ccc;
		  ddi = ddd;
		}
	      }
	    }
	  }
	}
	/* save result */
	rr[i] = rri;
	cc[i] = cci;
      }
    }
  }
}
