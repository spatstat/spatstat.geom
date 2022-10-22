/*

  fardist.h

  Code template for fardist.c

  Macros used:
      FNAME  function name
      SQUARED  #defined if squared distances should be returned.

  Copyright (C) Adrian Baddeley, Rolf Turner and Ege Rubak 2014
  Licence: GPL >= 2

  $Revision: 1.4 $  $Date: 2022/10/21 10:43:01 $


*/

void FNAME(
  /* inputs */
  int *nx,
  double *x0,
  double *xstep,   /* pixel grid dimensions */
  int *ny,
  double *y0,
  double *ystep,    /* pixel grid dimensions */
  int *np,
  double *xp,
  double *yp,      /* data points */
  /* outputs */
  double *dfar     /* output grid */
) { 
  int Nxcol, Nyrow, Npoints;
  int i, j, k, ijpos;
  double  X0, Y0, Xstep, Ystep, yi, xj;
  double d2, d2max, dx, dy;

  Nxcol   = *nx;
  Nyrow   = *ny;
  Npoints = *np;
  X0      = *x0;
  Y0      = *y0;
  Xstep   = *xstep;
  Ystep   = *ystep;

  if(Npoints == 0)
    return;

  /* loop over pixels */

  for(j = 0, xj = X0; j < Nxcol; j++, xj += Xstep) {

    R_CheckUserInterrupt();
    
    for(i = 0, yi = Y0; i < Nyrow; i++, yi += Ystep) {

      d2max = 0.0;
      
      for(k = 0; k < Npoints; k++) {
	
	dx = xj - xp[k];
	dy = yi - yp[k]; 
	d2 = dx * dx + dy * dy;
	if(d2 > d2max) 
	  d2max = d2;
	
      }

      ijpos = i + j * Nyrow;

#ifdef SQUARED
      dfar[ijpos] = d2max;
#else
      dfar[ijpos] = sqrt(d2max);
#endif

    /* end of loop over grid points (i, j) */
    }
  }
}



