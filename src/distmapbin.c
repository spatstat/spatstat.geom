/*
       distmapbin.c

       Distance transform of a discrete binary image
       (8-connected or 24-connected path metric)
       
       $Revision: 1.12 $ $Date: 2023/08/28 06:27:24 $
       
  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2023
  Licence: GNU Public Licence >= 2

*/

#include <math.h>
#include "raster.h"
#include <R_ext/Utils.h>

void dist_to_bdry(Raster *d);

void shape_raster(Raster *ras, void *data,
		  double xmin, double ymin, double xmax, double ymax,
		  int nrow, int ncol, int mrow, int mcol);

/* define core algorithms, using template 'distmapbin.h' */

#define FNAME distmap_bin
#define CONNECT 8
#include "distmapbin.h"
#undef FNAME
#undef CONNECT

#define FNAME dist24map_bin
#define CONNECT 24
#include "distmapbin.h"
#undef FNAME
#undef CONNECT

/* R interface */

void distmapbin(
  int *connect,           /* connectivity: 8 or 24 */
  double *xmin,
  double *ymin,
  double *xmax,
  double *ymax,  	  /* x, y dimensions */
  int *nr,
  int *nc,                /* raster dimensions
	                   EXCLUDING margin of 1 on each side */
  int   *inp,              /* input:  binary image */
  double *distances,	/* output: distance to nearest point */
  double *boundary       /* output: distance to boundary of rectangle */
	/* all images must have identical dimensions including a margin of 1 on each side */
) {
	Raster data, dist, bdist;

	void distmap_bin(Raster *in, Raster *dist);
	void dist24map_bin(Raster *in, Raster *dist);

	shape_raster( &data, (void *) inp, *xmin,*ymin,*xmax,*ymax,
			    *nr+2, *nc+2, 1, 1);
	shape_raster( &dist, (void *) distances,*xmin,*ymin,*xmax,*ymax,
			   *nr+2,*nc+2,1,1);
	shape_raster( &bdist, (void *) boundary, *xmin,*ymin,*xmax,*ymax,
			   *nr+2,*nc+2,1,1);

	if(*connect != 24) {
	  distmap_bin(&data, &dist);
	} else {
	  dist24map_bin(&data, &dist);
	}

	dist_to_bdry(&bdist);
}	
