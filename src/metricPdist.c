/*
  metricPdist.c

  Distance transform of binary pixel image
  using arbitrary metric

  This code #includes metricPdist.h multiple times.

  $Revision: 1.5 $  $Date: 2022/10/22 09:29:51 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

 */

/* Once-only declarations: */

#include <math.h>
#include "raster.h"

void dist_to_bdry(Raster *d);
void shape_raster(Raster *ras, void *data,
		  double xmin, double ymin, double xmax, double ymax,
		  int nrow, int ncol, int mrow, int mcol);

#define UNDEFINED -1
#define Is_Defined(I) (I >= 0)
#define Is_Undefined(I) (I < 0)

/*

  Definitions for each metric

  For each definition we need the following macros:
  FNAME          Name of the function (that will be called from R)
  MARGLIST       List of typed arguments to FNAME specifying the metric
  MTEMPDECLARE   Declaration and initialisation of variables for use by metric
  METRIC         Expression for calculating the metric (x1,y1,x2,y2)
*/

/* 
   (1) Rectangular metric

   Unit ball is a rectangle with width 1 unit, height 'aspect' units.

   mdt  = metric distance transform
   P    = pixel image input
   O    = orthogonally oriented to axis
   rect = rectangular 
*/

#define FNAME mdtPOrect
#define MARGLIST double *aspect
#define MTEMPDECLARE double asp; asp=*aspect
#define METRIC(X,Y,XX,YY) rectdist(X,Y,XX,YY,asp)

double rectdist(double x, double y, double xx, double yy, double asp)
{
  double dx, dy, d;
  dx = x-xx;
  dy = (y-yy)/asp;
  if(dx < 0) dx = -dx;
  if(dy < 0) dy = -dy;
  d = (dx > dy)? dx :  dy;
  return d;
}

#include "metricPdist.h"

#undef FNAME
#undef MARGLIST
#undef MARGDECLARE
#undef MTEMPDECLARE 
#undef METRIC

/* 
   (2) Convex ball metric

   Unit ball is a symmetric convex polygon

   mdt  = metric distance transform
   P    = pixel image input
   conv = convex set
*/

#define FNAME mdtPconv
#define MARGLIST int *ns, double *sx, double *sy
#define MTEMPDECLARE int Ns; Ns=*ns; 
#define METRIC(X,Y,XX,YY) convdist(X,Y,XX,YY,Ns,sx,sy)

double convdist(
  double x,
  double y,
  double xx,
  double yy,
  int Ns,
  double *sx,
  double *sy
) {
  int k;
  double dx, dy, dk, d;
  dx = x-xx;
  dy = y-yy;
  d = 0.0;
  for(k = 0; k < Ns; k++) {
    dk = dx * sx[k] + dy * sy[k];
    if(dk > d) d = dk;
  }
  return d;
}

#include "metricPdist.h"

#undef FNAME
#undef MARGLIST
#undef MARGDECLARE
#undef MTEMPDECLARE 
#undef METRIC


