/*

  knndistance.c

  K-th Nearest Neighbour Distances between points

  Copyright (C) Adrian Baddeley, Jens Oehlschlaegel and Rolf Turner 2000-2026
  Licence: GNU Public Licence >= 2

  $Revision: 1.13 $     $Date: 2026/03/16 07:10:13 $

  Function definitions are #included from knndist.h and knnXdist.h

  THE FOLLOWING FUNCTIONS ASSUME THAT y IS SORTED IN ASCENDING ORDER 

  SINGLE LIST:
  knndsort     k-th nearest neighbour distances
  knnd2sort    k-th nearest neighbour distances^2
  knnwhich     k-th nearest neighbours
  knnsort      k-th nearest neighbours and their distances
  knnsort2     k-th nearest neighbours and their distances^2

  ONE LIST TO ANOTHER LIST:
  knnXdist     Nearest neighbour distance from one list to another
  knnXdist2    Nearest neighbour distance^2 from one list to another
  knnXwhich    Nearest neighbour ID from one list to another
  knnX         Nearest neighbour ID & distance from one list to another
  knnX2        Nearest neighbour ID & distance^2 from one list to another

  ONE LIST TO ANOTHER OVERLAPPING LIST:
  knnXEdist    Nearest neighbour distance from one list to another, overlapping
  knnXEdist2   Nearest neighbour distance^2 from one list to another, overlapp
  knnXEwhich   Nearest neighbour ID from one list to another, overlapping
  knnXE        Nearest neighbour ID & distance 
  knnXE2       Nearest neighbour ID & distance^2

  Functions beginning 'knnX' can be accessed using 'knnXinterface'

*/

#undef SPATSTAT_DEBUG

#include <R.h>
#include <R_ext/Utils.h>
#include <math.h>

#include "yesno.h"

double sqrt(double x);

/* macros used - initially undefined */

#undef FNAME
#undef DIST
#undef WHICH
#undef EXCLUDE
#undef SQUARED


/* THE FOLLOWING CODE ASSUMES THAT y IS SORTED IN ASCENDING ORDER */

/* ------------------- one point pattern X --------------------- */

/* 
   knndsort 

   nearest neighbours 1:kmax

   returns distances only

*/

#define FNAME knndsort
#define DIST
#include "knndist.h"
#undef DIST
#undef FNAME

/* 
   knnd2sort 

   nearest neighbours 1:kmax

   returns squared distances only

*/

#define FNAME knnd2sort
#define DIST
#define SQUARED
#include "knndist.h"
#undef SQUARED
#undef DIST
#undef FNAME

/* 
   knnwhich

   nearest neighbours 1:kmax

   returns identifiers only

*/

#define FNAME knnwhich
#define WHICH
#include "knndist.h"
#undef WHICH
#undef FNAME

/* 
   knnsort 

   nearest neighbours 1:kmax

   returns distances and indices

*/

#define FNAME knnsort
#define DIST
#define WHICH
#include "knndist.h"
#undef WHICH
#undef DIST
#undef FNAME

/* 
   knnsort2 

   nearest neighbours 1:kmax

   returns squared distances and indices

*/

#define FNAME knnsort2
#define DIST
#define WHICH
#define SQUARED
#include "knndist.h"
#undef SQUARED
#undef WHICH
#undef DIST
#undef FNAME

/* --------------- two distinct point patterns X and Y --------------- */

/* Turn off the debugging tracer in knnXdist.h */
#undef TRACER

/* 
   knnXdist

   returns distances only

*/

#define FNAME knnXdist
#define DIST
#include "knnXdist.h"
#undef DIST
#undef FNAME

/* 
   knnXdist2

   returns squared distances only

*/

#define FNAME knnXdist2
#define DIST
#define SQUARED
#include "knnXdist.h"
#undef SQUARED
#undef DIST
#undef FNAME

/* 
   knnXwhich

   returns identifiers only

*/

#define FNAME knnXwhich
#define WHICH
#include "knnXdist.h"
#undef WHICH
#undef FNAME

/* 
   knnX 

   returns distances and indices

*/

#define FNAME knnX
#define DIST
#define WHICH
#include "knnXdist.h"
#undef WHICH
#undef DIST
#undef FNAME

/* 
   knnX2 

   returns distances and indices

*/

#define FNAME knnX2
#define DIST
#define WHICH
#define SQUARED
#include "knnXdist.h"
#undef SQUARED
#undef WHICH
#undef DIST
#undef FNAME

/* --------------- overlapping point patterns X and Y --------------- */

/* 
   knnXEdist

   returns distances only

*/

#define FNAME knnXEdist
#define DIST
#define EXCLUDE
#include "knnXdist.h"
#undef EXCLUDE
#undef DIST
#undef FNAME

/* 
   knnXEdist2

   returns squared distances only

*/

#define FNAME knnXEdist2
#define DIST
#define EXCLUDE
#define SQUARED
#include "knnXdist.h"
#undef SQUARED
#undef EXCLUDE
#undef DIST
#undef FNAME

/* 
   knnXEwhich

   returns identifiers only

*/

#define FNAME knnXEwhich
#define WHICH
#define EXCLUDE
#include "knnXdist.h"
#undef EXCLUDE
#undef WHICH
#undef FNAME

/* 
   knnXE 

   returns distances and indices

*/

#define FNAME knnXE
#define DIST
#define WHICH
#define EXCLUDE
#include "knnXdist.h"
#undef EXCLUDE
#undef WHICH
#undef DIST
#undef FNAME

/* 
   knnXE2

   returns distances and indices

*/

#define FNAME knnXE2
#define DIST
#define WHICH
#define EXCLUDE
#define SQUARED
#include "knnXdist.h"
#undef SQUARED
#undef EXCLUDE
#undef WHICH
#undef DIST
#undef FNAME


/* >>>>>>>>> GENERAL INTERFACE <<<<<<<<<<<<<<<< */

/* general interface for two patterns */

void knnXinterface(
  /* inputs */
  int *n1, double *x1, double *y1, int *id1,
  int *n2, double *x2, double *y2, int *id2,
  int *kmax,
  /* options */
  int *exclude,
  int *wantdist,
  int *wantwhich,
  int *squared,
  /* outputs */
  double *nnd,
  int *nnwhich,
  /* input (upper bound) */
  double *huge
     /* some inputs + outputs are not used in all functions */
) {
  int ex, di, wh, sq;
  ex = (*exclude != 0);
  di = (*wantdist != 0);
  wh = (*wantwhich != 0);
  sq = (*squared != 0);
  if(!ex) {
    /* ------------ two distinct point patterns: id1, id2 ignored -------- */
    if(di && wh) {
      /* distance and index */
      if(sq) {
	/* squared distance and index */
	knnX2(n1, x1, y1, id1, n2, x2, y2, id2, kmax, nnd, nnwhich, huge);
      } else {
	/* distance and index */
	knnX(n1, x1, y1, id1, n2, x2, y2, id2, kmax, nnd, nnwhich, huge);
      }
    } else if(di) {
      /* distance only */
      if(sq) {
	/* squared distance */
	knnXdist2(n1, x1, y1, id1, n2, x2, y2, id2, kmax, nnd, nnwhich, huge);
      } else {
	/* distance */
	knnXdist(n1, x1, y1, id1, n2, x2, y2, id2, kmax, nnd, nnwhich, huge);
      }
    } else if(wh) {
      /* index only */
      knnXwhich(n1, x1, y1, id1, n2, x2, y2, id2, kmax, nnd, nnwhich, huge);
    } 
  } else {
    /* ------  two patterns with overlap: use id1, id2 to avoid -------- */
    if(di && wh) {
      /* distance and index */
      if(sq) {
	/* squared distance and index */
	knnXE2(n1, x1, y1, id1, n2, x2, y2, id2, kmax, nnd, nnwhich, huge);
      } else {
	/* distance and index */
	knnXE(n1, x1, y1, id1, n2, x2, y2, id2, kmax, nnd, nnwhich, huge);
      }
    } else if(di) {
      /* distance only */
      if(sq) {
	/* squared distance */
	knnXEdist2(n1, x1, y1, id1, n2, x2, y2, id2, kmax, nnd, nnwhich, huge);
      } else {
	/* distance */
	knnXEdist(n1, x1, y1, id1, n2, x2, y2, id2, kmax, nnd, nnwhich, huge);
      }
    } else if(wh) {
      /* index only */
      knnXEwhich(n1, x1, y1, id1, n2, x2, y2, id2, kmax, nnd, nnwhich, huge);
    } 
  }
}


