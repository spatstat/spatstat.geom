/*

  nngrid.c

  Nearest Neighbour Distances from a pixel grid to a point pattern

  Copyright (C) Adrian Baddeley, Jens Oehlschlaegel and Rolf Turner 2000-2013
  Licence: GNU Public Licence >= 2

  $Revision: 1.7 $     $Date: 2026/03/11 07:58:01 $

  Function body definition is #included from nngrid.h 

  THE FOLLOWING FUNCTIONS ASSUME THAT x IS SORTED IN ASCENDING ORDER 

*/

#undef SPATSTAT_DEBUG

#include <R.h>
#include <R_ext/Utils.h>
#include <math.h>

#include "yesno.h"

double sqrt(double x);

/* initialise macros */
#undef FNAME
#undef DIST
#undef WHICH
#undef SQUARED

/* THE FOLLOWING CODE ASSUMES THAT x IS SORTED IN ASCENDING ORDER */

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

/* 
   nnGw 

   returns indices only

*/

#define FNAME nnGw
#define WHICH
#include "nngrid.h"
#undef FNAME
#undef WHICH

/* 
   nnGd2w

   returns _squared_ distances and indices

*/

#define FNAME nnGd2w
#define DIST
#define WHICH
#define SQUARED
#include "nngrid.h"
#undef FNAME
#undef DIST
#undef WHICH
#undef SQUARED

/* 
   nnGd2

   returns _squared_ distances only

*/

#define FNAME nnGd2
#define DIST
#define SQUARED
#include "nngrid.h"
#undef FNAME
#undef DIST
#undef SQUARED

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
  int *squared,
  /* outputs */
  double *nnd,
  int *nnwhich,
  /* upper bound on pairwise distance */
  double *huge
) {
  int di, wh, sq;
  di = (*wantdist != 0);
  wh = (*wantwhich != 0);
  sq = (*squared != 0);
  if(di && wh) {
    if(sq) {
      /* squared distance and index */
      nnGd2w(nx, x0, xstep, ny, y0, ystep, np, xp, yp, nnd, nnwhich, huge);
    } else {
      /* distance and index */
      nnGdw(nx, x0, xstep, ny, y0, ystep, np, xp, yp, nnd, nnwhich, huge);
    }
  } else if(di) {
    if(sq) {
      /* squared distance only */
      nnGd2(nx, x0, xstep, ny, y0, ystep, np, xp, yp, nnd, nnwhich, huge);
    } else {
      /* distance only */
      nnGd(nx, x0, xstep, ny, y0, ystep, np, xp, yp, nnd, nnwhich, huge);
    }
  } else if(wh) {
    /* index only */
    nnGw(nx, x0, xstep, ny, y0, ystep, np, xp, yp, nnd, nnwhich, huge);
  }
}

