/*
  scan.c

  Scan transform

  $Revision: 1.6 $ $Date: 2022/10/22 09:29:51 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/
#include <R.h>
#include <math.h>
#include "raster.h"

void shape_raster(Raster *ras, void *data,
		  double xmin, double ymin, double xmax, double ymax,
		  int nrow, int ncol, int mrow, int mcol);

void
Cscantrans(
  double *x,
  double *y,	/* data points */
  int	npt,
  double  R,    /* radius */
  Raster *out	/* scan image */
) {
  int	i,j,k,l,m;
  double  d2, R2;
  int   rmin, rmax, cmin, cmax, Rrow, Rcol, lmin, lmax, mmin, mmax;

  /* initialise raster */
  Clear(*out,int,0);
  
  /* If the list of data points is empty, ... exit now */
  if(npt == 0) 
    return;

  R2 = R * R;
  cmin = out->cmin;
  cmax = out->cmax;
  rmin = out->rmin;
  rmax = out->rmax;

  /* disc size in rows/columns */
  Rrow = (int) ceil(R/(out->ystep));
  Rcol = (int) ceil(R/(out->xstep));
  if(Rrow < 1) Rrow = 1; 
  if(Rcol < 1) Rcol = 1;
	
  /* run through points */
  for(i = 0; i < npt; i++) {
    j = RowIndex(*out,y[i]);
    k = ColIndex(*out,x[i]);
    lmin = j - Rrow;  if(lmin < rmin) lmin = rmin; 
    lmax = j + Rrow;  if(lmax > rmax) lmax = rmax;
    mmin = k - Rcol;  if(mmin < cmin) mmin = cmin; 
    mmax = k + Rcol;  if(mmax > cmax) mmax = cmax;

    for(l = lmin; l <= lmax; l++) {
      for(m = mmin; m <= mmax; m++) {
	d2 = DistanceToSquared(x[i],y[i],*out,l,m);
	if(d2 <= R2) 
	  Entry(*out,l,m,int) += 1;
      }
    }
  }
}

/* R interface */

void scantrans(
  double *x,
  double *y,		/* input data points */
  int	*n,
  double *xmin,
  double *ymin,
  double *xmax,
  double *ymax,  	/* guaranteed bounding box */
  int *nr,
  int *nc,		/* desired raster dimensions */
  double *R,            /* radius */
  /* output array */
  int *counts	        /* number of R-close points */
) {
  Raster out;
  int nrow, ncol, npoints;
  double r;

  nrow = *nr;
  ncol = *nc;
  npoints = *n;
  r = *R;

  shape_raster( &out, (void *) counts,
		*xmin,*ymin,*xmax,*ymax,
		nrow, ncol, 0, 0);
  Cscantrans(x, y, npoints, r, &out);
}	
