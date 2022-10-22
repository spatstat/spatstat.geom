/*
  
  raster.c

  shape_raster()     initialise a Raster structure
  
  $Revision: 1.2 $ $Date: 2022/10/22 02:32:10 $

*/

#include <math.h>
#include "raster.h"

void 
shape_raster(
  /* the raster structure to be initialised */
  Raster *ras,      
 /* pointer to data storage for pixel values */
  void *data,
  /* range of GRID COORDS excluding margin */  
  double xmin,
  double ymin,
  double xmax,
  double ymax,
  /* absolute dimensions of storage array */  
  int nrow,
  int ncol,
  /* margins for working */
  int mrow,
  int mcol  
) {
	ras->data	= data;
	ras->nrow 	= nrow;
	ras->ncol 	= ncol;
	ras->length 	= nrow * ncol;
	ras->rmin	= mrow;
	ras->rmax	= nrow - mrow - 1;
	ras->cmin	= mcol;
	ras->cmax	= ncol - mcol - 1;
	ras->x0		= 
	ras->xmin	= xmin;
	ras->x1 	=
	ras->xmax	= xmax;
	ras->y0		=
	ras->ymin	= ymin;
	ras->y1		=
	ras->ymax	= ymax;
	ras->xstep	= (xmax-xmin)/(ncol - 2 * mcol - 1);
	ras->ystep	= (ymax-ymin)/(nrow - 2 * mrow - 1);
	/* Rprintf("xstep,ystep = %lf,%lf\n", ras->xstep,ras->ystep);  */
}
