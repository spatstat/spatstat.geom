/*
       connectgraph.c

       Connected component transform for a finite graph

       cocoGraph: connected component labels for a discrete graph
                   specified by a list of edges
       
       $Revision: 1.1 $ $Date: 2023/07/18 03:28:58 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2023
  Licence: GNU Public Licence >= 2

       
*/

#include <math.h>
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>

#include "yesno.h"

void cocoGraph(
     /* inputs */
  int *nv,         /* number of graph vertices */
  int *ne,         /* number of edges */
  int *ie,
  int *je,         /* vectors of indices of ends of each edge */ 
  /* output */
  int *label,      /* vector of component labels for each vertex */
                   /* Component label is lowest serial number of
		        any vertex in the connected component */
  int *status          /* 0 if OK, 1 if overflow */
) {
  int Nv, Ne, i, j, k, niter, labi, labj, changed;
  
  Nv = *nv;
  Ne = *ne;

  /* initialise labels */
  for(k = 0; k < Nv; k++)
    label[k] = k;

  for(niter = 0; niter < Nv; niter++) {
    R_CheckUserInterrupt();
    changed = NO;
    for(k = 0; k < Ne; k++) {
      i = ie[k];
      j = je[k];
      labi = label[i];
      labj = label[j];
      if(labi < labj) {
	label[j] = labi;
	changed = YES;
      } else if(labj < labi) {
	label[i] = labj;
	changed = YES;
      } 
    }
    if(!changed) {
      /* algorithm has converged */
      *status = 0;
      return;
    }
  }
  /* error exit */   
  *status = 1;
  return;
}
