/*
       exactdist.c

       Exact distance transform of a point pattern
       (used to estimate the empty space function F)
       
       $Revision: 1.20 $ $Date: 2022/10/22 09:29:51 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

       Author: Adrian Baddeley

       Sketch of functionality:
            the 'data' are a finite list of points in R^2 
	    (x,y coordinates) and the 'output' is a real valued 
	    image whose entries are distances, with the value for
	    each pixel equalling the distance from that pixel
	    to the nearest point of the data pattern.
       
       Routines:

            exact_dt_R()       interface to R
	    exact_dt()         implementation of distance transform
	    dist_to_bdry()     compute distance to edge of image frame
	    shape_raster()     initialise a Raster structure
                          
       The appropriate calling sequence for exact_dt_R() 
       is exemplified in 'exactdt.R'
     
*/
#undef DEBUG

#include <math.h>
#include "raster.h"

#ifdef DEBUG
#include <R.h>
#endif

void shape_raster(Raster *ras, void *data,
		  double xmin, double ymin, double xmax, double ymax,
		  int nrow, int ncol, int mrow, int mcol);

void
exact_dt(
  double *x,
  double *y,		/* data points */
  int	npt,
  Raster *dist,		/* exact distance to nearest point */
  Raster *index		/* which point x[i],y[i] is closest */
) {
	int	i,j,k,l,m;
	double	d;
	int	ii;
	double	dd;
	/*	double  bdiag; */
	
	    /* initialise rasters */
#define UNDEFINED -1
#define Is_Defined(I) (I >= 0)
#define Is_Undefined(I) (I < 0)
	
	Clear(*index,int,UNDEFINED)
		
	d = 2.0 * DistanceSquared(dist->xmin,dist->ymin,dist->xmax,dist->ymax); 
	Clear(*dist,double,d)

	  /* If the list of data points is empty, ... exit now */
	if(npt == 0) 
	  return;

	for(i = 0; i < npt; i++) {
		/* Rprintf("%ld -> (%lf,%lf)\n", i, x[i], y[i]); */
		j = RowIndex(*dist,y[i]);
		k = ColIndex(*dist,x[i]);
		/* if(!Inside(*dist,j,k))
			Rprintf("(%ld,%ld) out of bounds\n",j,k);
		else if (!Inside(*dist,j+1,k+1))
			Rprintf("(%ld+1,%ld+1) out of bounds\n",j,k);
		*/
		for(l = j; l <= j+1; l++) 
		for(m = k; m <= k+1; m++) {
			d = DistanceToSquared(x[i],y[i],*index,l,m);
			if(   Is_Undefined(Entry(*index,l,m,int))
			   || Entry(*dist,l,m,double) > d)
			{
				/* Rprintf("writing (%ld,%ld) -> %ld\t%lf\n", l,m,i,d); */
				Entry(*index,l,m,int) = i;
				Entry(*dist,l,m,double) = d;
				/* Rprintf("checking: %ld, %lf\n",
				       Entry(*index,l,m,int),
				       Entry(*dist,l,m,double));
				 */
			}
		}
	}
/*
	for(j = 0; j <= index->nrow; j++)
		for(k = 0; k <= index->ncol; k++)
			Rprintf("[%ld,%ld] %ld\t%lf\n",
			       j,k,Entry(*index,j,k,int),Entry(*dist,j,k,double));
*/			
	/* how to update the distance values */
	
#define COMPARE(ROW,COL,RR,CC) \
	d = Entry(*dist,ROW,COL,double); \
	ii = Entry(*index,RR,CC,int); \
	/* Rprintf(" %lf\t (%ld,%ld) |-> %ld\n", d, RR, CC, ii); */ \
	if(Is_Defined(ii) /* && ii < npt */ \
	   && Entry(*dist,RR,CC,double) < d) { \
	     dd = DistanceSquared(x[ii],y[ii],Xpos(*index,COL),Ypos(*index,ROW)); \
	     if(dd < d) { \
		/* Rprintf("(%ld,%ld) <- %ld\n", ROW, COL, ii); */ \
		Entry(*index,ROW,COL,int) = ii; \
		Entry(*dist,ROW,COL,double) = dd; \
		/* Rprintf("checking: %ld, %lf\n", Entry(*index,ROW,COL,int), Entry(*dist,ROW,COL,double)); */\
	     } \
	}


	/* bound on diagonal step distance */
	/*	bdiag = sqrt(index->xstep * index->xstep + index->ystep * index->ystep); */
	
	/* forward pass */

	for(j = index->rmin; j <= index->rmax; j++)
	for(k = index->cmin; k <= index->cmax; k++) {
		/* Rprintf("Neighbourhood of (%ld,%ld):\n", j,k); */
		COMPARE(j,k, j-1,k-1)
		COMPARE(j,k, j-1,  k)
		COMPARE(j,k, j-1,k+1)
		COMPARE(j,k, j,  k-1)
	}

	/* backward pass */

	for(j = index->rmax; j >= index->rmin; j--)
	for(k = index->cmax; k >= index->cmin; k--) {
		COMPARE(j,k, j+1,k+1)
		COMPARE(j,k, j+1,  k)
		COMPARE(j,k, j+1,k-1)
		COMPARE(j,k, j,  k+1)
	}

	/* take square roots of the distances^2 */

	for(j = index->rmin; j <= index->rmax; j++)
	for(k = index->cmin; k <= index->cmax; k++) 
	        Entry(*dist,j,k,double) = sqrt(Entry(*dist,j,k,double));
	
}	

#define MIN(A,B) (((A) < (B)) ? (A) : (B))

/* distance to frame boundary from each raster point */

void
dist_to_bdry(Raster *d)
{
  int j, k;
  double x, y, xd, yd, Xmin, Xmax, Ymin, Ymax;
  /* Frame limits */
  Xmin = d->xmin - d->xstep/2.0;
  Xmax = d->xmax + d->xstep/2.0;
  Ymin = d->ymin - d->ystep/2.0;
  Ymax = d->ymax + d->ystep/2.0;

#ifdef DEBUG
  Rprintf("xmin=%lf,xmax=%lf\nymin=%lf,ymax=%lf\n",
	  d->xmin, d->xmax, d->ymin, d->ymax);
  Rprintf("xstep=%lf,ystep=%lf\n",
	  d->xstep, d->ystep);
  Rprintf("Xmin=%lf,Xmax=%lf\nYmin=%lf,Ymax=%lf\n",
	  Xmin,Xmax, Ymin, Ymax);
  Rprintf("(Xpos,Ypos)(cmin,rmin) = (%lf, %lf)\n",
	  Xpos(*d,cmin), Ypos(*d,rmin));
#endif
  
  for(j = d->rmin; j <= d->rmax;j++) {
    y = Ypos(*d,j);
    yd = MIN(y - Ymin, Ymax - y);
    for(k = d->cmin; k <= d->cmax;k++) {
      x = Xpos(*d,k);
      xd = MIN(x - Xmin, Xmax - x);
      Entry(*d,j,k,double) = MIN(xd,yd);
    }
  }
}

/* R interface */

void exact_dt_R(
  double *x,
  double *y,		/* input data points */
  int	*npt,
  double *xmin,
  double *ymin,
  double *xmax,
  double *ymax,  	/* guaranteed bounding box */
  int *nr,
  int *nc,		/* desired raster dimensions
			   EXCLUDING margins */
  int *mr,
  int *mc,           /* margins */
  /* output arrays */
  double *distances,	/* distance to nearest point */
  int   *indices,	        /* index to nearest point */
  double *boundary	/* distance to boundary */
) {
	Raster dist, index, bdist;
	int mrow, mcol, nrow, ncol;

	mrow = *mr;
	mcol = *mc;

	/* full dimensions */
	nrow = *nr + 2 * mrow;
	ncol = *nc + 2 * mcol;
	
	shape_raster( &dist, (void *) distances,*xmin,*ymin,*xmax,*ymax,
		      nrow, ncol, mrow, mcol);
	shape_raster( &index, (void *) indices, *xmin,*ymin,*xmax,*ymax,
		      nrow, ncol, mrow, mcol);
	shape_raster( &bdist, (void *) boundary, *xmin,*ymin,*xmax,*ymax,
		      nrow, ncol, mrow, mcol);
	
	exact_dt(x, y, (int) *npt, &dist, &index);
	dist_to_bdry(&bdist);
}	
