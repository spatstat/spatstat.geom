/*

  knngrid.c

  K-th Nearest Neighbour Distances from a pixel grid to a point pattern

  Copyright (C) Adrian Baddeley, Jens Oehlschlaegel and Rolf Turner 2000-2022
  Licence: GNU Public Licence >= 2

  $Revision: 1.8 $     $Date: 2022/10/22 02:44:26 $

  Function body definition is #included from knngrid.h 

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
   knnGdw

   nearest neighbours 1:kmax

   returns distances and indices

*/

#define FNAME knnGdw
#define DIST
#define WHICH
#include "knngrid.h"
#undef FNAME
#undef DIST
#undef WHICH

/* 
   knnGd

   nearest neighbours 1:kmax

   returns distances only

*/

#define FNAME knnGd
#define DIST
#include "knngrid.h"
#undef FNAME
#undef DIST
#undef WHICH

/* 
   knnGw 

   nearest neighbours 1:kmax

   returns indices only

*/

#define FNAME knnGw
#define WHICH
#include "knngrid.h"
#undef FNAME
#undef DIST
#undef WHICH

/* >>>>>>>>>>> GENERAL INTERFACE <<<<<<<<<<<<<<<< */

/* general interface */

void knnGinterface(
  /* inputs */
  int *nx, double *x0, double *xstep,
  int *ny, double *y0, double *ystep,  /* pixel grid dimensions */
  int *np, double *xp, double *yp,     /* data points */
  int *kmax,
  /* options */
  int *wantdist,
  int *wantwhich,
  /* outputs */
  double *nnd,
  int *nnwhich,
  /* upper bound on pairwise distance */
  double *huge
  /* some inputs + outputs are not used in all functions */
) {
  int di, wh;
  di = (*wantdist != 0);
  wh = (*wantwhich != 0);
  if(di && wh) {
    knnGdw(nx, x0, xstep, ny, y0, ystep, np, xp, yp, kmax, nnd, nnwhich, huge);
  } else if(di) {
    knnGd(nx, x0, xstep, ny, y0, ystep, np, xp, yp, kmax, nnd, nnwhich, huge);
  } else if(wh) {
    knnGw(nx, x0, xstep, ny, y0, ystep, np, xp, yp, kmax, nnd, nnwhich, huge);
  }
}

