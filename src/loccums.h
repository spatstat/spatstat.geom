/*
  loccums.h

  C template for loccum.c

  data-to-data functions

  $Revision: 1.7 $ $Date: 2022/10/21 10:43:01 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

  macros: 

  FNAME    function name
  NULVAL   initial value (empty sum = 0, empty product = 1)
  INC(A,B) increment operation A += B or A *= B

*/

void FNAME(
  /* inputs */
  int *n,
  double *x,
  double *y,
  double *v,
  int *nr,
  double *rmax, 
  /* output */
  double *ans  /* matrix of column vectors of functions for each point */
) {
  int N, Nr, Nans;
  double Rmax;

  int i, j, k, kmin, maxchunk, columnstart;
  double Rmax2, rstep, xi, yi;
  double dx, dy, dx2, d2, d, contrib;

  N    = *n;
  Nr   = *nr;
  Rmax = *rmax;

  if(N == 0) 
    return;

  rstep = Rmax/(Nr-1);
  Rmax2 = Rmax * Rmax;
  Nans  = Nr * N;

  /* initialise products to 1 */
  OUTERCHUNKLOOP(k, Nans, maxchunk, 8196) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(k, Nans, maxchunk, 8196) {
      ans[k] = NULVAL;
    }
  }
   
  OUTERCHUNKLOOP(i, N, maxchunk, 8196) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, N, maxchunk, 8196) {
      xi = x[i];
      yi = y[i];
      columnstart = Nr * i; /* start position for f_i(.) in 'ans' */
      /* 
	 process backward until |dx| > Rmax
      */
      if(i > 0) {
	for(j=i-1; j >= 0; j--) {
	  dx = x[j] - xi;
	  dx2 = dx * dx;
	  if(dx2 > Rmax2) 
	    break;
	  dy = y[j] - yi;
	  d2 = dx2 + dy * dy;
	  if(d2 <= Rmax2) {
	    d = sqrt(d2);
	    kmin = (int) ceil(d/rstep);
	    if(kmin < Nr) {
	      contrib = v[j];
	      for(k = kmin; k < Nr; k++) 
		INC(ans[columnstart + k] , contrib);
	    }
	  }
	}
      }
      /* 
	 process forward until |dx| > Rmax
      */
      if(i < N - 1) {
	for(j=i+1; j < N; j++) {
	  dx = x[j] - xi;
	  dx2 = dx * dx;
	  if(dx2 > Rmax2) 
	    break;
	  dy = y[j] - yi;
	  d2 = dx2 + dy * dy;
	  if(d2 <= Rmax2) {
	    d = sqrt(d2);
	    kmin = (int) ceil(d/rstep);
	    if(kmin < Nr) {
	      contrib = v[j];
	      for(k = kmin; k < Nr; k++) 
		INC(ans[columnstart + k] , contrib);
	    }
	  }
	}
      }
    }   
  }
}
