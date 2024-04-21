/*

  bdrymask.c

  Boundary pixels of binary mask

  Copyright (C) Adrian Baddeley, Rolf Turner and Ege Rubak 2014
  Licence: GPL >= 2

  $Revision: 1.4 $  $Date: 2022/10/20 10:57:43 $


*/

#include <R.h>
#include <R_ext/Utils.h>
#include <math.h>

void bdrymask(
  /* inputs */
  int *nx,
  int *ny,
  int *m,
  /* outputs */
  int *b
) { 
  int Nxcol, Nyrow, Nx1, Ny1;
  int i, j, mij;

  Nxcol   = *nx;
  Nyrow   = *ny;
  Nx1 = Nxcol - 1;
  Ny1 = Nyrow - 1;

#define MAT(A,I,J) A[(I) + (J) * Nyrow]

  /* loop over pixels */

  for(j = 0; j < Nxcol; j++) {

    R_CheckUserInterrupt();
    
    for(i = 0; i < Nyrow; i++) {

      mij = MAT(m, i, j);
      if(i == 0 || i == Ny1 || j == 0 || j == Nx1) {
	MAT(b, i, j) = mij;
      } else if((mij != MAT(m, (i-1), j)) ||
		(mij != MAT(m, (i+1), j)) ||
		(mij != MAT(m, i, (j-1))) ||
		(mij != MAT(m, i, (j+1)))) {
	MAT(b, i, j) = 1;
      }
    }
  }
}



