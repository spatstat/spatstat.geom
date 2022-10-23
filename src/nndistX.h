
#if (1 == 0)
/*
  nndistX.h

  Code template for C functions supporting nncross

  THE FOLLOWING CODE ASSUMES THAT LISTS ARE SORTED
  IN ASCENDING ORDER OF y COORDINATE

  This code is #included multiple times in nndistance.c 
  Variables used:
        FNAME     function name
        DIST      #defined if function returns distance to nearest neighbour
	WHICH     #defined if function returns id of nearest neighbour
        EXCLUDE   #defined if exclusion mechanism is used
  Either or both DIST and WHICH may be defined.

  When EXCLUDE is defined,
  code numbers id1, id2 are attached to the patterns X and Y respectively, 
  such that
  x1[i], y1[i] and x2[j], y2[j] are the same point iff id1[i] = id2[j].

  Copyright (C) Adrian Baddeley, Jens Oehlschlagel and Rolf Turner 2000-2012
  Licence: GPL >= 2

  $Revision: 1.8 $  $Date: 2022/10/23 05:21:52 $


*/
#endif

#undef USEJ

#if ((defined WHICH) | (!defined EXCLUDE))
#define USEJ
#endif

void FNAME(
  int *n1, double *x1, double *y1, int *id1, 
  int *n2, double *x2, double *y2, int *id2, 
  double *nnd,
  int *nnwhich, 
  double *huge
  /* some inputs + outputs are not used in all functions */
) { 
  int npoints1, npoints2, maxchunk, i, jleft, jright, lastjwhich;
  double d2, d2min, x1i, y1i, dx, dy, dy2, hu, hu2;
#ifdef USEJ
  int jwhich;
#endif  
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
  
  /* loop in chunks of 2^16 */

  i = 0; maxchunk = 0; 
  while(i < npoints1) {

    R_CheckUserInterrupt();

    maxchunk += 65536; 
    if(maxchunk > npoints1) maxchunk = npoints1;

    for(; i < maxchunk; i++) {

      d2min = hu2;
      x1i = x1[i];
      y1i = y1[i];
#ifdef EXCLUDE
      id1i = id1[i];
#endif
#ifdef USEJ
      jwhich = -1;
#endif

      if(lastjwhich < npoints2) {  /* always true if EXCLUDE is defined */
	/* search forward from previous nearest neighbour  */
	for(jright = lastjwhich; jright < npoints2; ++jright)
	  {
	    dy = y2[jright] - y1i;
	    dy2 = dy * dy; 
	    if(dy2 > d2min) /* note that dy2 >= d2min could break too early */
	      break;
#ifdef EXCLUDE
	    /* do not compare identical points */
	    if(id2[jright] != id1i) {
#endif
	      dx = x2[jright] - x1i;
	      d2 =  dx * dx + dy2;
	      if (d2 < d2min) {
		d2min = d2;
#ifdef USEJ	
		jwhich = jright;
#endif		
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
	    dy = y1i - y2[jleft];
	    dy2 = dy * dy;
	    if(dy2 > d2min) /* note that dy2 >= d2min could break too early */
	      break;
#ifdef EXCLUDE
	    /* do not compare identical points */
	    if(id2[jleft] != id1i) {
#endif
	      dx = x2[jleft] - x1i;
	      d2 =  dx * dx + dy2;
	      if (d2 < d2min) {
		d2min = d2;
#ifdef USEJ	
		jwhich = jleft;
#endif		
	      }
#ifdef EXCLUDE
	    }
#endif
	  }
	/* end backward search */
      }
      /* commit values */
#ifdef DIST
      nnd[i] = sqrt(d2min);
#endif
#ifdef WHICH
      nnwhich[i] = jwhich + 1; /* R indexing */
#endif

#ifndef EXCLUDE
      lastjwhich = jwhich;
#endif
      
    }
  }
}
