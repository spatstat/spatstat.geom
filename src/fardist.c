/*

  fardist.c

  Furthest data point from each grid point

  Uses code template 'fardist.h'

  Copyright (C) Adrian Baddeley, Rolf Turner and Ege Rubak 2014
  Licence: GPL >= 2

  $Revision: 1.3 $  $Date: 2022/10/22 09:29:51 $


*/

#include <R.h>
#include <R_ext/Utils.h>
#include <math.h>

double sqrt(double x);

#define FNAME fardistgrid
#undef  SQUARED
#include "fardist.h"
#undef FNAME

#define FNAME fardist2grid
#define  SQUARED
#include "fardist.h"
