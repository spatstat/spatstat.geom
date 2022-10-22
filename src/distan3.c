/*

  distan3.c

  Distances between pairs of 3D points

  $Revision: 1.7 $     $Date: 2022/10/22 02:44:42 $

  D3pairdist      Pairwise distances
  D3pair2dist     Pairwise distances squared
  D3pairPdist     Pairwise distances with periodic correction
  D3pairP2dist    Pairwise distances squared, with periodic correction

  D3crossdist     Pairwise distances for two sets of points
  D3cross2dist    Pairwise distances squared, for two sets of points
  D3crossPdist    Pairwise distances for two sets of points, periodic correction

  matchxyz       Find matches between two sets of points   

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

 */

#include <math.h>
/* #include <stdio.h> */

double sqrt();

void D3pairdist(
  /* inputs */
  int *n,
  double *x,
  double *y,
  double *z,
  int *squared,
  /* output */
  double *d
) { 
  void D3pair1dist(int *n, double *x, double *y, double *z, double *d);
  void D3pair2dist(int *n, double *x, double *y, double *z, double *d);
  if(*squared == 0) {
    D3pair1dist(n, x, y, z, d);
  } else {
    D3pair2dist(n, x, y, z, d);
  }
}


void D3pair1dist(
  /* inputs */
  int *n,
  double *x,
  double *y,
  double *z,
  /* output */
  double *d
) { 
  int i, j, npoints; 
  double *dp;
  double xi, yi, zi, dx, dy, dz, dist;

  npoints = *n;

  /* set d[0,0] = 0 */
  *d = 0.0;

  for (i=1; i < npoints; i++) 
    {
      xi = x[i];
      yi = y[i];
      zi = z[i];
      /* point at the start of column i */
      dp = d + i * npoints;
      /* set diagonal to zero */
      dp[i] = 0.0;
      for (j=0; j < i; j++)
	{
	  dx = x[j] - xi;
	  dy = y[j] - yi;
	  dz = z[j] - zi;
	  dist = sqrt( dx * dx + dy * dy + dz * dz ); 
	  /* upper triangle */
	  *dp = dist;
	  ++dp;
	  /* lower triangle */
	  d[ j * npoints + i] = dist;
	}
    }
}

/* squared distances */

void D3pair2dist(
  /* inputs */
  int *n,
  double *x,
  double *y,
  double *z,
  /* output */
  double *d
) { 
  int i, j, npoints; 
  double *dp;
  double xi, yi, zi, dx, dy, dz, dist;

  npoints = *n;

  /* set d[0,0] = 0 */
  *d = 0.0;

  for (i=1; i < npoints; i++) 
    {
      xi = x[i];
      yi = y[i];
      zi = z[i];
      /* point at the start of column i */
      dp = d + i * npoints;
      /* set diagonal to zero */
      dp[i] = 0.0;
      for (j=0; j < i; j++)
	{
	  dx = x[j] - xi;
	  dy = y[j] - yi;
	  dz = z[j] - zi;
	  dist = dx * dx + dy * dy + dz * dz; 
	  /* upper triangle */
	  *dp = dist;
	  ++dp;
	  /* lower triangle */
	  d[ j * npoints + i] = dist;
	}
    }
}

void D3crossdist(
  /* inputs */
  int *nfrom, double *xfrom, double *yfrom, double *zfrom,
  int *nto,   double *xto, double *yto, double *zto,
  int *squared,
  /* output */
  double *d
) {
  void D3cross1dist(int *nfrom, double *xfrom, double *yfrom, double *zfrom,
		    int *nto,   double *xto, double *yto, double *zto,
		    double *d);
  void D3cross2dist(int *nfrom, double *xfrom, double *yfrom, double *zfrom,
		    int *nto,   double *xto, double *yto, double *zto,
		    double *d);
  if(*squared == 0) {
    D3cross1dist(nfrom, xfrom, yfrom, zfrom, nto, xto, yto, zto, d);
  } else {
    D3cross2dist(nfrom, xfrom, yfrom, zfrom, nto, xto, yto, zto, d);
  }
}


void D3cross1dist(
  /* inputs */
  int *nfrom, double *xfrom, double *yfrom, double *zfrom,
  int *nto,   double *xto, double *yto, double *zto,
  /* output */
  double *d
) {
  int i, j, nf, nt; 
  double *dptr;
  double xj, yj, zj, dx, dy, dz;

  nf = *nfrom;
  nt = *nto;

  dptr = d;

  for (j=0; j < nt; j++) {
    xj = xto[j];
    yj = yto[j];
    zj = zto[j];
    for(i = 0; i < nf; i++, dptr++) {
	dx = xj - xfrom[i];
	dy = yj - yfrom[i];
	dz = zj - zfrom[i];
	*dptr = sqrt( dx * dx + dy * dy + dz * dz ); 
    }
  }
}

/* squared distances */

void D3cross2dist(
  /* inputs */
  int *nfrom, double *xfrom, double *yfrom, double *zfrom,
  int *nto,   double *xto, double *yto, double *zto,
  /* output */
  double *d
) {
  int i, j, nf, nt; 
  double *dptr;
  double xj, yj, zj, dx, dy, dz;

  nf = *nfrom;
  nt = *nto;

  dptr = d;

  for (j=0; j < nt; j++) {
    xj = xto[j];
    yj = yto[j];
    zj = zto[j];
    for(i = 0; i < nf; i++, dptr++) {
	dx = xj - xfrom[i];
	dy = yj - yfrom[i];
	dz = zj - zfrom[i];
	*dptr = dx * dx + dy * dy + dz * dz; 
    }
  }
}



/* distances with periodic correction */

void D3pairPdist(
  /* inputs */
  int *n,
  double *x,
  double *y,
  double *z,
  double *xwidth,
  double *yheight,
  double *zdepth,
  int *squared,
  /* output */
  double *d
) {
  void D3pairP1dist(int *n, double *x, double *y, double *z,
		    double *xwidth, double *yheight, double *zdepth,
		    double *d);
  void D3pairP2dist(int *n, double *x, double *y, double *z,
		    double *xwidth, double *yheight, double *zdepth,
		    double *d);
  if(*squared == 0) {
    D3pairP1dist(n, x, y, z, xwidth, yheight, zdepth, d);
  } else {
    D3pairP2dist(n, x, y, z, xwidth, yheight, zdepth, d);
  }
}

void D3pairP1dist(
  /* inputs */
  int *n,
  double *x,
  double *y,
  double *z,
  double *xwidth,
  double *yheight,
  double *zdepth,
  /* output */
  double *d
) {
  int i, j, npoints; 
  double *dp;
  double xi, yi, zi, dx, dy, dz, dx2, dy2, dz2, dx2p, dy2p, dz2p, dist, wide, high, deep;

  npoints = *n;
  wide = *xwidth;
  high = *yheight;
  deep = *zdepth;

  /* set d[0,0] = 0 */
  *d = 0.0;

  for (i=1; i < npoints; i++) 
    {
      xi = x[i];
      yi = y[i];
      zi = z[i];
      /* point at the start of column i */
      dp = d + i * npoints;
      /* set diagonal to zero */
      dp[i] = 0.0;
      for (j=0; j < i; j++)
	{
	  dx = x[j] - xi;
	  dy = y[j] - yi;
	  dz = z[j] - zi;
	  dx2p = dx * dx;
	  dy2p = dy * dy;
	  dz2p = dz * dz;
	  dx2 = (dx - wide) * (dx - wide);
	  dy2 = (dy - high) * (dy - high);
	  dz2 = (dz - deep) * (dz - deep);
	  if(dx2 < dx2p) dx2p = dx2;
	  if(dy2 < dy2p) dy2p = dy2;
	  if(dz2 < dz2p) dz2p = dz2;
	  dx2 = (dx + wide) * (dx + wide);
	  dy2 = (dy + high) * (dy + high);
	  dz2 = (dz + deep) * (dz + deep);
	  if(dx2 < dx2p) dx2p = dx2;
	  if(dy2 < dy2p) dy2p = dy2;
	  if(dz2 < dz2p) dz2p = dz2;
	  dist = sqrt( dx2p + dy2p + dz2p ); 
	  /* upper triangle */
	  *dp = dist;
	  ++dp;
	  /* lower triangle */
	  d[ j * npoints + i] = dist;
	}
    }
}

/* same function without the sqrt */

void D3pairP2dist(
  /* inputs */
  int *n,
  double *x,
  double *y,
  double *z,
  double *xwidth,
  double *yheight,
  double *zdepth,
  /* output */
  double *d
) {
  int i, j, npoints; 
  double *dp;
  double xi, yi, zi, dx, dy, dz, dx2, dy2, dz2, dx2p, dy2p, dz2p, dist, wide, high, deep;

  npoints = *n;
  wide = *xwidth;
  high = *yheight;
  deep = *zdepth;

  /* set d[0,0] = 0 */
  *d = 0.0;

  for (i=1; i < npoints; i++) 
    {
      xi = x[i];
      yi = y[i];
      zi = z[i];
      /* point at the start of column i */
      dp = d + i * npoints;
      /* set diagonal to zero */
      dp[i] = 0.0;
      for (j=0; j < i; j++)
	{
	  dx = x[j] - xi;
	  dy = y[j] - yi;
	  dz = z[j] - zi;
	  dx2p = dx * dx;
	  dy2p = dy * dy;
	  dz2p = dz * dz;
	  dx2 = (dx - wide) * (dx - wide);
	  dy2 = (dy - high) * (dy - high);
	  dz2 = (dz - deep) * (dz - deep);
	  if(dx2 < dx2p) dx2p = dx2;
	  if(dy2 < dy2p) dy2p = dy2;
	  if(dz2 < dz2p) dz2p = dz2;
	  dx2 = (dx + wide) * (dx + wide);
	  dy2 = (dy + high) * (dy + high);
	  dz2 = (dz + deep) * (dz + deep);
	  if(dx2 < dx2p) dx2p = dx2;
	  if(dy2 < dy2p) dy2p = dy2;
	  if(dz2 < dz2p) dz2p = dz2;
	  dist = dx2p + dy2p + dz2p; 
	  /* upper triangle */
	  *dp = dist;
	  ++dp;
	  /* lower triangle */
	  d[ j * npoints + i] = dist;
	}
    }
}

void D3crossPdist(
  /* inputs */
  int *nfrom, double *xfrom, double *yfrom, double *zfrom,
  int *nto, double *xto, double *yto, double *zto,
  double *xwidth, double *yheight, double *zdepth,
  int *squared,
  /* output */
  double *d
) {
  void D3crossP1dist(int *nfrom, double *xfrom, double *yfrom, double *zfrom,
		     int *nto, double *xto, double *yto, double *zto,
		     double *xwidth, double *yheight, double *zdepth,
		     double *d);
  void D3crossP2dist(int *nfrom, double *xfrom, double *yfrom, double *zfrom,
		     int *nto, double *xto, double *yto, double *zto,
		     double *xwidth, double *yheight, double *zdepth,
		     double *d);
  if(*squared == 0) {
    D3crossP1dist(nfrom, xfrom, yfrom, zfrom, 
		  nto, xto, yto, zto, 
		  xwidth, yheight, zdepth, 
		  d);
  } else {
    D3crossP2dist(nfrom, xfrom, yfrom, zfrom, 
		  nto, xto, yto, zto, 
		  xwidth, yheight, zdepth, 
		  d);
  }
}


void D3crossP1dist(
  /* inputs */
  int *nfrom, double *xfrom, double *yfrom, double *zfrom,
  int *nto, double *xto, double *yto, double *zto,
  double *xwidth, double *yheight, double *zdepth,
  /* output */
  double *d
) { 
  int i, j, nf, nt; 
  double *dptr;
  double xj, yj, zj, dx, dy, dz, dx2, dy2, dz2, dx2p, dy2p, dz2p, wide, high, deep;

  nf = *nfrom;
  nt = *nto;
  wide = *xwidth;
  high = *yheight;
  deep = *zdepth;

  dptr = d;

  for (j=0; j < nt; j++) {
    xj = xto[j];
    yj = yto[j];
    zj = zto[j];
    for(i = 0; i < nf; i++, dptr++) {
	dx = xj - xfrom[i];
	dy = yj - yfrom[i];
	dz = zj - zfrom[i];
	dx2p = dx * dx;
	dy2p = dy * dy;
	dz2p = dz * dz;
	dx2 = (dx - wide) * (dx - wide);
	dy2 = (dy - high) * (dy - high);
	dz2 = (dz - deep) * (dz - deep);
	if(dx2 < dx2p) dx2p = dx2;
	if(dy2 < dy2p) dy2p = dy2;
	if(dz2 < dz2p) dz2p = dz2;
	dx2 = (dx + wide) * (dx + wide);
	dy2 = (dy + high) * (dy + high);
	dz2 = (dz + deep) * (dz + deep);
	if(dx2 < dx2p) dx2p = dx2;
	if(dy2 < dy2p) dy2p = dy2;
	if(dz2 < dz2p) dz2p = dz2;
	*dptr = sqrt( dx2p + dy2p + dz2p ); 
    }
  }
}


void D3crossP2dist(
  /* inputs */
  int *nfrom, double *xfrom, double *yfrom, double *zfrom,
  int *nto, double *xto, double *yto, double *zto,
  double *xwidth, double *yheight, double *zdepth,
  /* output */
  double *d
) { 
  int i, j, nf, nt; 
  double *dptr;
  double xj, yj, zj, dx, dy, dz, dx2, dy2, dz2, dx2p, dy2p, dz2p, wide, high, deep;

  nf = *nfrom;
  nt = *nto;
  wide = *xwidth;
  high = *yheight;
  deep = *zdepth;

  dptr = d;

  for (j=0; j < nt; j++) {
    xj = xto[j];
    yj = yto[j];
    zj = zto[j];
    for(i = 0; i < nf; i++, dptr++) {
	dx = xj - xfrom[i];
	dy = yj - yfrom[i];
	dz = zj - zfrom[i];
	dx2p = dx * dx;
	dy2p = dy * dy;
	dz2p = dz * dz;
	dx2 = (dx - wide) * (dx - wide);
	dy2 = (dy - high) * (dy - high);
	dz2 = (dz - deep) * (dz - deep);
	if(dx2 < dx2p) dx2p = dx2;
	if(dy2 < dy2p) dy2p = dy2;
	if(dz2 < dz2p) dz2p = dz2;
	dx2 = (dx + wide) * (dx + wide);
	dy2 = (dy + high) * (dy + high);
	dz2 = (dz + deep) * (dz + deep);
	if(dx2 < dx2p) dx2p = dx2;
	if(dy2 < dy2p) dy2p = dy2;
	if(dz2 < dz2p) dz2p = dz2;
	*dptr = dx2p + dy2p + dz2p; 
    }
  }
}

/*

  matchxyz

  Find matches between two lists of points

 */

void matchxyz(
  /* inputs */
  int *na,  double *xa, double *ya, double *za,
  int *nb,  double *xb, double *yb, double *zb,
  /* output */
  int *match
) { 
  int i, j, Na, Nb; 
  double xai, yai, zai;

  Na = *na;
  Nb = *nb;

  for (i=1; i < Na; i++) 
    {
      xai = xa[i];
      yai = ya[i];
      zai = za[i];
      match[i] = 0;
      for (j=0; j < Nb; j++) 
	if(xai == xb[j] && yai == yb[j] && zai == zb[i]) {
	  match[i] = j;
	  break;
	}
    }
}

