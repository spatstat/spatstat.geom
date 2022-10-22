/*

  nn3DdistX.h

  Code template for nearest-neighbour algorithms for 3D point patterns

  Input is two point patterns - supports 'nncross'

  This code is #included multiple times in nn3Ddist.c
  Variables used:
        FNAME     function name
        DIST      #defined if function returns distance to nearest neighbour
	WHICH     #defined if function returns id of nearest neighbour
	EXCLUDE   #defined if the two patterns may include common points
	          (which are not to be counted as neighbours)

  Either or both DIST and WHICH may be defined.

  THE FOLLOWING CODE ASSUMES THAT BOTH POINT PATTERNS ARE SORTED
  IN ASCENDING ORDER OF THE z COORDINATE

  If EXCLUDE is #defined, 
   Code numbers id1, id2 are attached to the patterns X and Y respectively, 
   such that
   x1[i], y1[i] and x2[j], y2[j] are the same point iff id1[i] = id2[j].

  $Revision: 1.8 $ $Date: 2022/10/22 09:22:51 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/

void FNAME(
  /* first point pattern */	   
  int *n1,
  double *x1,
  double *y1,
  double *z1,
  int *id1, 
  /* second point pattern */	   
  int *n2,
  double *x2,
  double *y2,
  double *z2,
  int *id2,
  /* outputs */
  double *nnd,     /* n.n. distances */
  int *nnwhich,    /* n.d. identifiers */
  /* prior upper bound on pairwise distances */
  double *huge
) { 
  int npoints1, npoints2, i, j, jwhich, lastjwhich;
  double d2, d2min, x1i, y1i, z1i, dx, dy, dz, dz2, hu, hu2;
#ifdef EXCLUDE
  int id1i;
#endif

  hu = *huge;
  hu2 = hu * hu;

  npoints1 = *n1;
  npoints2 = *n2;

  if(npoints1 == 0 || npoints2 == 0)
    return;

  lastjwhich = 0;  /* remains unchanged if EXCLUDE is defined */

  for(i = 0; i < npoints1; i++) {
    
    R_CheckUserInterrupt();
    
    d2min = hu2;
    jwhich = -1;
    x1i = x1[i];
    y1i = y1[i];
    z1i = z1[i];
#ifdef EXCLUDE
    id1i = id1[i];
#endif

    /* search backward from previous nearest neighbour */
    if(lastjwhich > 0) {  /* always true if EXCLUDE is defined */
      for(j = lastjwhich - 1; j >= 0; --j) {
	dz = z2[j] - z1i;
	dz2 = dz * dz;
	if(dz2 > d2min)
	  break;
#ifdef EXCLUDE
	/* do not compare identical points */
	if(id2[j] != id1i) {
#endif
	  dx = x2[j] - x1i;
	  dy = y2[j] - y1i;
	  d2 =  dx * dx + dy * dy + dz2;
	  if (d2 < d2min) {
	    d2min = d2;
	    jwhich = j;
	  }
#ifdef EXCLUDE
	}
#endif
      }
    }

    /* search forward from previous nearest neighbour  */
    if(lastjwhich < npoints2) { /* always false if EXCLUDE is defined */
      for(j = lastjwhich; j < npoints2; ++j) {
	dz = z2[j] - z1i;
	dz2 = dz * dz;
	if(dz2 > d2min)
	  break;
#ifdef EXCLUDE
	/* do not compare identical points */
	if(id2[j] != id1i) {
#endif
	  dx = x2[j] - x1i;
	  dy = y2[j] - y1i;
	  d2 =  dx * dx + dy * dy + dz2;
	  if (d2 < d2min) {
	    d2min = d2;
	    jwhich = j;
	  }
#ifdef EXCLUDE
	}
#endif
      }
    }
#ifdef DIST
    nnd[i] = sqrt(d2min);
#endif
#ifdef WHICH
    /* convert to R indexing */
    nnwhich[i] = jwhich + 1;
#endif
#ifndef EXCLUDE    
    lastjwhich = jwhich;
#endif    
  }
}
