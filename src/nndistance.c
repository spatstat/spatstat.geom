/*

  nndistance.c

  Nearest Neighbour Distances between points

  Copyright (C) Adrian Baddeley, Jens Oehlschlaegel and Rolf Turner 2000-2026
  Licence: GNU Public Licence >= 2

  $Revision: 1.27 $     $Date: 2026/03/16 02:29:37 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2024
  Licence: GNU Public Licence >= 2

  THE FOLLOWING FUNCTIONS ASSUME THAT y IS SORTED IN ASCENDING ORDER 

  SINGLE LIST:
  nndistsort    Nearest neighbour distances 
  nndist2sort   Nearest neighbour distances^2
  nnwhichsort   Nearest neighbours
  nnsort        Nearest neighbours & distances
  nnsort2        Nearest neighbours & distances^2

  ONE LIST TO ANOTHER LIST:
  nnXdist       Nearest neighbour distance from one list to another
  nnXdist2      Nearest neighbour distance^2 from one list to another
  nnXwhich      Nearest neighbour ID from one list to another
  nnX           Nearest neighbour ID & distance from one list to another
  nnX2          Nearest neighbour ID & distance^2 from one list to another

  ONE LIST TO ANOTHER OVERLAPPING LIST:
  nnXEdist      Nearest neighbour distance from one list to another, overlapping
  nnXEdist2     Nearest neighbour distance^2 from one list to another, overlap
  nnXEwhich     Nearest neighbour ID from one list to another, overlapping
  nnXE          Nearest neighbour ID & distance 
  nnXE2          Nearest neighbour ID & distance^2

  Functions beginning 'nnX' can be accessed through 'nnXinterface'
*/

#undef SPATSTAT_DEBUG

#include <R.h>
#include <R_ext/Utils.h>
#include <math.h>

#include "yesno.h"

double sqrt(double x);

/* 

  This file #includes the files 'nndist.h' and 'nndistX.h' multiple times
  with different values of the following macros.

*/

#undef FNAME
#undef DIST
#undef WHICH
#undef SQUARED
#undef EXCLUDE


/* THE FOLLOWING CODE ASSUMES THAT y IS SORTED IN ASCENDING ORDER */


/* ------------------- one point pattern X --------------------- */

/* 
   nndistsort: nearest neighbour distances 
*/

#define FNAME nndistsort
#define DIST
#include "nndist.h"
#undef DIST
#undef FNAME

/* 
   nndist2sort: squared nearest neighbour distances
*/

#define FNAME nndist2sort
#define DIST
#define SQUARED
#include "nndist.h"
#undef SQUARED
#undef DIST
#undef FNAME

/* 
   nnwhichsort: id of nearest neighbour 
*/

#define FNAME nnwhichsort
#define WHICH
#include "nndist.h"
#undef WHICH
#undef FNAME

/* 
   nnsort: distance & id of nearest neighbour 
*/

#define FNAME nnsort
#define DIST
#define WHICH
#include "nndist.h"
#undef WHICH
#undef DIST
#undef FNAME

/* 
   nnsort2: squared distance & id of nearest neighbour 
*/

#define FNAME nnsort2
#define DIST
#define WHICH
#define SQUARED
#include "nndist.h"
#undef SQUARED
#undef WHICH
#undef DIST
#undef FNAME

/* --------------- two distinct point patterns X and Y  ----------------- */


/* 
   nnXdist:  nearest neighbour distance
	      (from each point of X to the nearest point of Y)
*/

#define FNAME nnXdist
#define DIST
#include "nndistX.h"
#undef DIST
#undef FNAME

/* 
   nnXdist2:  squared nearest neighbour distance
	      (from each point of X to the nearest point of Y)
*/

#define FNAME nnXdist2
#define DIST
#define SQUARED
#include "nndistX.h"
#undef SQUARED
#undef DIST
#undef FNAME

/* 
   nnXwhich:  nearest neighbour id
*/

#define FNAME nnXwhich
#define WHICH
#include "nndistX.h"
#undef WHICH
#undef FNAME

/* 
   nnX:  nearest neighbour distance and id
*/

#define FNAME nnX
#define DIST
#define WHICH
#include "nndistX.h"
#undef WHICH
#undef DIST
#undef FNAME

/* 
   nnX2:  squared nearest neighbour distance and id
*/

#define FNAME nnX2
#define DIST
#define WHICH
#define SQUARED
#include "nndistX.h"
#undef SQUARED
#undef WHICH
#undef DIST
#undef FNAME

/* --------------- two point patterns X and Y with common points --------- */

/*
   Code numbers id1, id2 are attached to the patterns X and Y respectively, 
   such that
   x1[i], y1[i] and x2[j], y2[j] are the same point iff id1[i] = id2[j].
*/

/* 
   nnXEdist:  similar to nnXdist
          but allows X and Y to include common points
          (which are not to be counted as neighbours)
*/

#define FNAME nnXEdist
#define DIST
#define EXCLUDE
#include "nndistX.h"
#undef EXCLUDE
#undef DIST
#undef FNAME

/* 
   nnXEdist2: similar to nnXdist2 (squared distances)
          but allows X and Y to include common points
          (which are not to be counted as neighbours)
*/

#define FNAME nnXEdist2
#define DIST
#define EXCLUDE
#define SQUARED
#include "nndistX.h"
#undef SQUARED
#undef EXCLUDE
#undef DIST
#undef FNAME

/* 
   nnXEwhich:  similar to nnXwhich
          but allows X and Y to include common points
          (which are not to be counted as neighbours)
*/

#define FNAME nnXEwhich
#define WHICH
#define EXCLUDE
#include "nndistX.h"
#undef EXCLUDE
#undef WHICH
#undef FNAME

/* 
   nnXE:  similar to nnX
          but allows X and Y to include common points
          (which are not to be counted as neighbours)
*/

#define FNAME nnXE
#define DIST
#define WHICH
#define EXCLUDE
#include "nndistX.h"
#undef EXCLUDE
#undef WHICH
#undef DIST
#undef FNAME

/* 
   nnXE2:  similar to nnX2 (squared distances)
          but allows X and Y to include common points
          (which are not to be counted as neighbours)
*/

#define FNAME nnXE2
#define DIST
#define WHICH
#define EXCLUDE
#define SQUARED
#include "nndistX.h"
#undef SQUARED
#undef EXCLUDE
#undef WHICH
#undef DIST
#undef FNAME


/* general interface */

void nnXinterface(
   /* first point pattern */
   int *n1,
   double *x1,
   double *y1,
   int *id1, 
   /* second point pattern */
   int *n2,
   double *x2,
   double *y2,
   int *id2, 
   /* options */
   int *exclude,
   int *wantdist,
   int *wantwhich,
   int *squared,
   /* outputs */
   double *nnd,
   int *nnwhich,
   /* largest possible distance */
   double *huge
) {
  /* defined above */
  /* void nnX(), nnXdist(), nnXwhich(); */
  /* void nnXE(), nnXEdist(), nnXEwhich(); */
  int ex, di, wh, sq;
  ex = (*exclude != 0);
  di = (*wantdist != 0);
  wh = (*wantwhich != 0);
  sq = (*squared != 0);
  if(!ex) {
    /* two distinct point patterns: id1, id2 ignored */ 
    if(di && wh) {
      /* both distance and index */
      if(sq) {
	/* squared distance and index */
	nnX2(n1, x1, y1, id1, n2, x2, y2, id2, nnd, nnwhich, huge);
      } else {
	/* distance and index */
	nnX(n1, x1, y1, id1, n2, x2, y2, id2, nnd, nnwhich, huge);
      }
    } else if(di) {
      /* distance only */
      if(sq) {
	/* squared distance */
	nnXdist2(n1, x1, y1, id1, n2, x2, y2, id2, nnd, nnwhich, huge);
      } else {
	/* distance */
	nnXdist(n1, x1, y1, id1, n2, x2, y2, id2, nnd, nnwhich, huge);
      }
    } else if(wh) {
      /* index only */
      nnXwhich(n1, x1, y1, id1, n2, x2, y2, id2, nnd, nnwhich, huge);
    } 
  } else {
    /* overlapping patterns: use id1, id2 to avoid overlap */
    if(di && wh) {
      /* both distance and index */
      if(sq) {
	/* squared distance and index */
	nnXE2(n1, x1, y1, id1, n2, x2, y2, id2, nnd, nnwhich, huge);
      } else {
	/* distance and index */
	nnXE(n1, x1, y1, id1, n2, x2, y2, id2, nnd, nnwhich, huge);
      }
    } else if(di) {
      /* distance only */
      if(sq) {
	/* squared distance */
	nnXEdist2(n1, x1, y1, id1, n2, x2, y2, id2, nnd, nnwhich, huge);
      } else {
	/* distance */
	nnXEdist(n1, x1, y1, id1, n2, x2, y2, id2, nnd, nnwhich, huge);
      }
    } else if(wh) {
      /* index only */
      nnXEwhich(n1, x1, y1, id1, n2, x2, y2, id2, nnd, nnwhich, huge);
    } 
  }
}

