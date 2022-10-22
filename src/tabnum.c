/*
  tabnum.c

  table(x) or tapply(x, w, sum)
  where x is numeric and we are given the sorted unique values

  $Revision: 1.5 $ $Date: 2022/10/22 02:32:10 $

*/

#include <math.h>
#include <R_ext/Utils.h>
#include "chunkloop.h"

void tabnum(
  int *nx, 
  double *x,  /* values (sorted) */
  int *nv, 
  double *v,  /* unique values (sorted) */
  double *z   /* output */
) { 
  int i, j, Nx, Nv, maxchunk; 
  double xi;
  
  Nx = *nx;
  Nv = *nv;

  j = 0;
  
  OUTERCHUNKLOOP(i, Nx, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, Nx, maxchunk, 16384) {
      xi = x[i];
      /* Find the smallest v[j] greater than or equal to x[i] */
      for( ; j < Nv && xi > v[j]; j++)
	;
      /* increment */
      if(j < Nv)
	z[j] += 1.0;
    }
  }
}

void tabsumweight(
     int *nx, 
     double *x,  /* values */
     double *w,  /* weights */
     int *nv, 
     double *v,  /* unique values (sorted) */
     double *z   /* output */
) { 
  int i, j, Nx, Nv, maxchunk; 
  double xi;
  
  Nx = *nx;
  Nv = *nv;

  j = 0;
  
  OUTERCHUNKLOOP(i, Nx, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(i, Nx, maxchunk, 16384) {
      xi = x[i];
      /* Find the smallest v[j] greater than or equal to x[i] */
      for(; j < Nv && xi > v[j]; j++)
	;
      /* add weight */
      if(j < Nv)
	z[j] += w[i];
    }
  }
}

