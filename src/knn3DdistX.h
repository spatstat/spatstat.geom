
#if (1 == 0)
/*
  knn3DdistX.h

  Code template for C functions supporting nncross 
  for k-nearest neighbours (k > 1)
  for 3D point patterns

  THE FOLLOWING CODE ASSUMES THAT LISTS ARE SORTED
  IN ASCENDING ORDER OF z COORDINATE

  This code is #included multiple times in nn3Ddist.c 
  Variables used:
        FNAME     function name
        DIST      #defined if function returns distance to nearest neighbour
	WHICH     #defined if function returns id of nearest neighbour
        EXCLUDE   #defined if exclusion mechanism is used
  Either or both DIST and WHICH may be defined.

  When EXCLUDE is defined,
  code numbers id1, id2 are attached to the patterns X and Y respectively, 
  such that X[i] and Y[j] are the same point iff id1[i] = id2[j].

  Copyright (C) Adrian Baddeley, Jens Oehlschlagel and Rolf Turner 2000-2013
  Licence: GPL >= 2

  $Revision: 1.5 $  $Date: 2023/05/09 04:59:28 $


  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/
#endif

void FNAME(
  /* inputs */
  int *n1, double *x1, double *y1, double *z1, int *id1,
  int *n2, double *x2, double *y2, double *z2, int *id2,
  int *kmax,
  /* output matrices (n1 * kmax) in ROW MAJOR order */
  double *nnd,
  int *nnwhich,
  double *huge
  /* some inputs + outputs are not used in all functions */
) { 
  int npoints1, npoints2, nk, nk1;
  int maxchunk, i, jleft, jright, lastjwhich, unsorted, k, k1;
  double d2, d2minK, x1i, y1i, z1i, dx, dy, dz, dz2, hu, hu2, tmp;
  double *d2min; 
#ifdef WHICH
  int *which;
  int itmp;
#endif
#ifdef EXCLUDE
  int id1i;
#else
  int jwhich;
#endif

  npoints1 = *n1;
  npoints2 = *n2;
  nk       = *kmax;
  nk1      = nk - 1;
  hu       = *huge;
  hu2      = hu * hu;

  if(npoints1 == 0 || npoints2 == 0)
    return;

  lastjwhich = 0; /* remains unchanged if EXCLUDE is defined */

  /* 
     create space to store the nearest neighbour distances and indices
     for the current point
  */

  d2min = (double *) R_alloc((size_t) nk, sizeof(double));
#ifdef WHICH
  which = (int *) R_alloc((size_t) nk, sizeof(int));
#endif

  /* loop in chunks of 2^16 */

  i = 0; maxchunk = 0; 
  while(i < npoints1) {

    R_CheckUserInterrupt();

    maxchunk += 65536; 
    if(maxchunk > npoints1) maxchunk = npoints1;

    for(; i < maxchunk; i++) {

      /* initialise nn distances and indices */
      d2minK = hu2;
#ifndef EXCLUDE      
      jwhich = -1;
#endif      
      for(k = 0; k < nk; k++) {
	d2min[k] = hu2;
#ifdef WHICH
	which[k] = -1;
#endif
      }

      x1i = x1[i];
      y1i = y1[i];
      z1i = z1[i];
#ifdef EXCLUDE
      id1i = id1[i];
#endif

      if(lastjwhich < npoints2) { /* always true if EXCLUDE is defined */
	/* search forward from previous nearest neighbour  */
	for(jright = lastjwhich; jright < npoints2; ++jright)
	  {
	    dz = z2[jright] - z1i;
	    dz2 = dz * dz; 
	    if(dz2 > d2minK) /* note that dz2 >= d2minK could break too early */
	      break;
#ifdef EXCLUDE
	    /* do not compare identical points */
	    if(id2[jright] != id1i) {
#endif
	      dy = y2[jright] - y1i;
	      d2 =  dy * dy + dz2;
	      if(d2 < d2minK) {
		dx = x2[jright] - x1i;
		d2 =  dx * dx + d2;
		if (d2 < d2minK) {
		  /* overwrite last entry in list of neighbours */
		  d2min[nk1] = d2;
#ifndef EXCLUDE
		  jwhich = jright;
#endif		  
#ifdef WHICH
		  which[nk1] = jright;
#endif
		  /* bubble sort */
		  unsorted = YES;
		  for(k = nk1; unsorted && k > 0; k--) {
		    k1 = k - 1;
		    if(d2min[k] < d2min[k1]) {
		      /* swap entries */
		      tmp  = d2min[k1];
		      d2min[k1] = d2min[k];
		      d2min[k] = tmp;
#ifdef WHICH
		      itmp = which[k1];
		      which[k1] = which[k];
		      which[k] = itmp;
#endif
		    } else {
		      unsorted = NO;
		    }
		  }
		  /* adjust maximum distance */
		  d2minK = d2min[nk1];
		}
	      }
#ifdef EXCLUDE
	    }
#endif
	  }
	/* end forward search */
      }
      if(lastjwhich > 0) { /* always false if EXCLUDE is defined */
	/* search backward from previous nearest neighbour */
	for(jleft = lastjwhich - 1; jleft >= 0; --jleft)
	  {
	    dz = z1i - z2[jleft];
	    dz2 = dz * dz;
	    if(dz2 > d2minK) /* note that dz2 >= d2minK could break too early */
	      break;
#ifdef EXCLUDE
	    /* do not compare identical points */
	    if(id2[jleft] != id1i) {
#endif
	      dy = y2[jleft] - y1i;
	      d2 =  dy * dy + dz2;
	      if(d2 < d2minK) {
		dx = x2[jleft] - x1i;
		d2 =  dx * dx + d2;
		if (d2 < d2minK) {
		  /* overwrite last entry in list of neighbours */
		  d2min[nk1] = d2;
#ifndef EXCLUDE
		  jwhich = jleft;
#endif		  
#ifdef WHICH
		  which[nk1] = jleft;
#endif
		  /* bubble sort */
		  unsorted = YES;
		  for(k = nk1; unsorted && k > 0; k--) {
		    k1 = k - 1;
		    if(d2min[k] < d2min[k1]) {
		      /* swap entries */
		      tmp  = d2min[k1];
		      d2min[k1] = d2min[k];
		      d2min[k] = tmp;
#ifdef WHICH
		      itmp = which[k1];
		      which[k1] = which[k];
		      which[k] = itmp;
#endif
		    } else {
		      unsorted = NO;
		    }
		  }
		  /* adjust maximum distance */
		  d2minK = d2min[nk1];
		}
	      }
#ifdef EXCLUDE
	    }
#endif
	  }
	/* end backward search */
      }
      /* copy nn distances for point i 
	 to output matrix in ROW MAJOR order
      */
      for(k = 0; k < nk; k++) {
#ifdef DIST
	nnd[nk * i + k] = sqrt(d2min[k]);
#endif
#ifdef WHICH
	nnwhich[nk * i + k] = which[k] + 1;  /* R indexing */
#endif
      }
#ifndef EXCLUDE      
      /* save index of last neighbour encountered */
      lastjwhich = jwhich;
#endif      
      /* end of loop over points i */
    }
  }
}

