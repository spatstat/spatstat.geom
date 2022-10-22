/*

  nngrid.c

  Nearest Neighbour Distances from a pixel grid to a point pattern

  Copyright (C) Adrian Baddeley, Jens Oehlschlaegel and Rolf Turner 2000-2013
  Licence: GNU Public Licence >= 2

  $Revision: 1.5 $     $Date: 2022/10/22 02:32:10 $

  Function body definition is #included from nngrid.h 

  THE FOLLOWING FUNCTIONS ASSUME THAT x IS SORTED IN ASCENDING ORDER 

*/

#undef SPATSTAT_DEBUG

#include <R.h>
#include <R_ext/Utils.h>
#include <math.h>

#include "yesno.h"

double sqrt();

/* THE FOLLOWING CODE ASSUMES THAT x IS SORTED IN ASCENDING ORDER */

#undef FNAME
#undef DIST
#undef WHICH

/* 
   nnGdw

   returns distances and indices

*/

#define FNAME nnGdw
#define DIST
#define WHICH
#include "nngrid.h"
#undef FNAME
#undef DIST
#undef WHICH

/* 
   nnGd

   returns distances only

*/

#define FNAME nnGd
#define DIST
#include "nngrid.h"
#undef FNAME
#undef DIST
#undef WHICH

/* 
   nnGw 

   returns indices only

*/

#define FNAME nnGw
#define WHICH
#include "nngrid.h"
#undef FNAME
#undef DIST
#undef WHICH

/* general interface */

void nnGinterface(
  /* pixel grid dimensions */
  int *nx,
  double *x0,
  double *xstep,  
  int *ny,
  double *y0,
  double *ystep,
  /* data points */
  int *np,
  double *xp,
  double *yp,
  /* options */  
  int *wantdist,
  int *wantwhich, 
  /* outputs */
  double *nnd,
  int *nnwhich,
  /* upper bound on pairwise distance */
  double *huge
) {
  int di, wh;
  di = (*wantdist != 0);
  wh = (*wantwhich != 0);
  if(di && wh) {
    nnGdw(nx, x0, xstep, ny, y0, ystep, np, xp, yp, nnd, nnwhich, huge);
  } else if(di) {
    nnGd(nx, x0, xstep, ny, y0, ystep, np, xp, yp, nnd, nnwhich, huge);
  } else if(wh) {
    nnGw(nx, x0, xstep, ny, y0, ystep, np, xp, yp, nnd, nnwhich, huge);
  }
}

