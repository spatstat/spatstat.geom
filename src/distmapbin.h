/*
       distmapbin.h

       Distance transform of a discrete binary image

       Template for core algorithm

       This file is #included in 'distmapbin.c' several times
       with different values of the macros:

       FNAME    Name of function
       CONNECT  Connectivity (8 or 24)

       $Revision: 1.2 $ $Date: 2023/08/28 07:37:38 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2023
  Licence: GNU Public Licence >= 2

*/

#define DISTANCE(ROW, COL) Entry(*dist, ROW, COL, double)
#define MASKTRUE(ROW, COL) (Entry(*in, ROW, COL, int) != 0)
#define MASKFALSE(ROW, COL) (Entry(*in, ROW, COL, int) == 0)
#define UPDATE(D, ROW, COL, STEP) \
	dnew = STEP + DISTANCE(ROW, COL); \
        if(D > dnew) D = dnew

void
FNAME(
    Raster  *in,            /* input:  binary image */
    Raster  *dist           /* output: distance to nearest point */
	/* rasters must have been dimensioned by shape_raster()
	   and must all have identical dimensions and margins */
) {
	int	j,k;
	double	d, dnew;
	double  xstep, ystep, diagstep, huge;
	int rmin, rmax, cmin, cmax;
	
#if (CONNECT == 24)
	double  LstepXXy, LstepYYx;
#endif
	
	/* distances between neighbouring pixels */
	xstep = in->xstep;
	ystep = in->ystep;
	if(xstep < 0) xstep = -xstep;
	if(ystep < 0) ystep = -ystep;
	diagstep = sqrt(xstep * xstep + ystep * ystep);
#if (CONNECT == 24)
	LstepXXy = sqrt(4.0 * xstep * xstep + ystep * ystep);
	LstepYYx = sqrt(xstep * xstep + 4.0 * ystep * ystep);
#endif	
	/* effectively infinite distance */
	huge = 2.0 * Distance(dist->xmin,dist->ymin,dist->xmax,dist->ymax); 

	/* image boundaries */
	rmin = in->rmin;
	rmax = in->rmax;
	cmin = in->cmin;
	cmax = in->cmax;

	/* initialise edges to boundary condition */
	for(j = rmin-1; j <= rmax+1; j++) {
	  DISTANCE(j, cmin-1) = (MASKTRUE(j, cmin-1)) ? 0.0 : huge;
	  DISTANCE(j, cmax+1) = (MASKTRUE(j, cmax+1)) ? 0.0 : huge;
	}
	for(k = cmin-1; k <= cmax+1; k++) {
	  DISTANCE(rmin-1, k) = (MASKTRUE(rmin-1, k)) ? 0.0 : huge;
	  DISTANCE(rmax+1, k) = (MASKTRUE(rmax+1, k)) ? 0.0 : huge;
	}

	/* forward pass */

	for(j = rmin; j <= rmax; j++) {
	  R_CheckUserInterrupt();
	  for(k = cmin; k <= cmax; k++) {
	    if(MASKTRUE(j, k))
	      d = DISTANCE(j, k) = 0.0;
	    else {
	      d = huge;
	      UPDATE(d, j-1, k-1, diagstep);
	      UPDATE(d, j-1,   k, ystep);
	      UPDATE(d, j-1, k+1, diagstep);
	      UPDATE(d,   j, k-1, xstep);
#if (CONNECT == 24)
	      if(j > rmin) {
		UPDATE(d, j-2, k-1, LstepYYx);
		UPDATE(d, j-2, k+1, LstepYYx);
	      }
	      if(k > cmin) {
		UPDATE(d, j-1, k-2, LstepXXy);
	      }
	      if(k < cmax) {
		UPDATE(d, j-1, k+2, LstepXXy);
	      }
#endif
	      DISTANCE(j,k) = d;
	    }
	  }
	}

	/* backward pass */

	for(j = rmax; j >= rmin; j--) {
	  R_CheckUserInterrupt();
	  for(k = cmax; k >= cmin; k--) {
	    if(MASKFALSE(j,k)) {
	      d = DISTANCE(j,k);
	      UPDATE(d, j+1, k+1, diagstep);
	      UPDATE(d, j+1,   k, ystep);
	      UPDATE(d, j+1, k-1, diagstep);
	      UPDATE(d,   j, k+1, xstep);
#if (CONNECT == 24)
	      if(j < rmax) {
		UPDATE(d, j+2, k-1, LstepYYx);
		UPDATE(d, j+2, k+1, LstepYYx);
	      }
	      if(k > cmin) {
		UPDATE(d, j+1, k-2, LstepXXy);
	      }
	      if(k < cmax) {
		UPDATE(d, j+1, k+2, LstepXXy);
	      }
#endif
	      DISTANCE(j,k) = d;
	    } 
	  }
	}
}

#undef DISTANCE
#undef MASKTRUE
#undef MASKFALSE
#undef UPDATE

