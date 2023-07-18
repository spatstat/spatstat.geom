/* 
   connectpix.h

   Code template for connected component transform of image 

   This file is #included multiple times in 'connectpix.c' using macros
	   INAME  Internal workhorse function name
           RNAME  R interface function name
           STYPE  Pixel value storage type (int or double)
           CONN   connectivity (4 or 8)

   $Revision: 1.6 $ $Date: 2023/07/18 03:45:27 $

   Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2023
   Licence: GNU Public Licence >= 2

*/

/* First the internal workhorse function */

void
INAME(Raster  *im)
     /* raster must have been dimensioned by shape_raster() */
     /* Pixel values assumed to be 0 in background, and 
        distinct nonzero values in foreground */
{
  int	j,k;
  int rmin, rmax, cmin, cmax;
  STYPE label, curlabel, minlabel;
  int anychanged;

  /* image boundaries */
  rmin = im->rmin;
  rmax = im->rmax;
  cmin = im->cmin;
  cmax = im->cmax;

  anychanged = 1;

#define ENTRY(ROW, COL) Entry(*im, ROW, COL, STYPE)

#define UPDATE(ROW,COL,BEST,NEW) \
     NEW = ENTRY(ROW, COL); \
     if(NEW != 0 && NEW < BEST) \
       BEST = NEW

  while(anychanged != 0) {
    anychanged = 0;
    R_CheckUserInterrupt();
    for(j = rmin; j <= rmax; j++) {
      for(k = cmin; k <= cmax; k++) {
	curlabel = ENTRY(j, k);
	if(curlabel != 0) {
	  minlabel = curlabel;
#if (CONN == 8)	  
	  UPDATE(j-1, k-1, minlabel, label);
#endif	  
	  UPDATE(j-1, k,   minlabel, label);
#if (CONN == 8)	  
	  UPDATE(j-1, k+1, minlabel, label);
#endif	  
	  UPDATE(j,   k-1, minlabel, label);
	  UPDATE(j,   k,   minlabel, label);
	  UPDATE(j,   k+1, minlabel, label);
#if (CONN == 8)	  
	  UPDATE(j+1, k-1, minlabel, label);
#endif	  
	  UPDATE(j+1, k,   minlabel, label);
#if (CONN == 8)	  
	  UPDATE(j+1, k+1, minlabel, label);
#endif	  
	  if(minlabel < curlabel) {
	    ENTRY(j, k) = minlabel;
	    anychanged = 1;
	  }
	}
      }
    }
  }
}

/* The R interface function */

void RNAME(
  STYPE *mat,        /* input */
  int *nr,
  int *nc          /* raster dimensions
			   EXCLUDING margin of 1 on each side */
) {
  Raster im;

  shape_raster( &im, (void *) mat, 
		(double) 1, (double) 1,
		(double) *nc, (double) *nr, 
		*nr+2, *nc+2, 1, 1);
  INAME(&im);
}	


#undef ENTRY
#undef UPDATE

