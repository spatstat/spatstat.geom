/*
       connectpix.c

       Connected component transforms for pixel image

       coco8int  8-connected, integer image
       coco4int  4-connected, integer image
       coco8dbl  8-connected, double precision image
       coco4dbl  4-connected, double precision image

       The input is a pixel image in which background pixels have value 0
       and all non-background pixels have been initialised to different nonzero 
       values. 

       (Double precision images are needed when the raster is too large
       for every pixel to be labelled with a different integer.)

       This code repeatedly #includes 'connectpix.h' 
       
       $Revision: 1.18 $ $Date: 2023/07/18 04:04:05 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2023
  Licence: GNU Public Licence >= 2

       
*/

#include <math.h>
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>

#include "raster.h"
#include "yesno.h"

void shape_raster(Raster *ras, void *data,
		  double xmin, double ymin, double xmax, double ymax,
		  int nrow, int ncol, int mrow, int mcol);


/* 
   integer, 8-connected 
*/

#define STYPE int
#define CONN 8
#define INAME Iconcom8
#define RNAME coco8int
#include "connectpix.h"
#undef INAME
#undef RNAME
#undef CONN
#undef STYPE

/* 
   integer, 4-connected 
*/

#define STYPE int
#define CONN 4
#define INAME Iconcom4
#define RNAME coco4int
#include "connectpix.h"
#undef INAME
#undef RNAME
#undef CONN
#undef STYPE

/* 
   double, 8-connected 
*/

#define STYPE double
#define CONN 8
#define INAME Dconcom8
#define RNAME coco8dbl
#include "connectpix.h"
#undef INAME
#undef RNAME
#undef CONN
#undef STYPE

/* 
   double, 4-connected 
*/

#define STYPE double
#define CONN 4
#define INAME Dconcom4
#define RNAME coco4dbl
#include "connectpix.h"
#undef INAME
#undef RNAME
#undef CONN
#undef STYPE



