#include <R.h>
#include <R_ext/Utils.h>
#include "chunkloop.h"
#include "looptest.h"

/*

  Epent.c

  $Revision: 1.4 $     $Date: 2026/01/08 10:19:30 $

  C implementation of 'eval' for Penttinen interaction potential

  Assumes point patterns are sorted in increasing order of x coordinate

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2026
  Licence: GNU Public Licence >= 2

*/

double sqrt(double x);
double acos(double x);
double log(double x);

void Epent(
  /* inputs */
  /* query points */
  int *nnsource,  
  double *xsource,
  double *ysource,
  int *idsource,
  /* data points */
  int *nntarget,  
  double *xtarget,
  double *ytarget,
  int *idtarget,
  /* model parameters */
  double *radius,
  /* output */
  double *values
) {
  int nsource, ntarget, maxchunk, j, i, ileft, idsourcej;
  double xsourcej, ysourcej, xleft, dx, dy, dx2, d2;
  double rad, reach, reach2, reach2pluseps;
  double z, z2, omz2, contrib, total;

  nsource = *nnsource;
  ntarget = *nntarget;
  rad     = *radius;

  if(nsource == 0 || ntarget == 0) 
    return;

  reach         = 2.0 * rad;
  reach2        = reach * reach;
  reach2pluseps = reach2 + EPSILON(reach2);

  ileft = 0;

  OUTERCHUNKLOOP(j, nsource, maxchunk, 65536) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(j, nsource, maxchunk, 65536) {
      total = 0.0;
      xsourcej = xsource[j];
      ysourcej = ysource[j];
      idsourcej = idsource[j];

      /* 
	 adjust starting point
      */

      xleft  = xsourcej - reach;
      while((xtarget[ileft] < xleft) && (ileft+1 < ntarget))
	++ileft;

      /* 
	 process until dx > reach 
      */
      for(i=ileft; i < ntarget; i++) {
	dx = xtarget[i] - xsourcej;
	dx2 = dx * dx;
	if(dx2 > reach2pluseps) 
	  break;
	if(idtarget[i] != idsourcej) {
	  /* valid pair */
	  dy = ytarget[i] - ysourcej;
	  d2 = dx2 + dy * dy;
	  if(d2 <= reach2) {
	    /* close pair */
	    z2 = d2/reach2;
	    omz2 = 1.0 - z2;
	    if(omz2 >= 0.0) {
	      z = sqrt(z2);
	      contrib = acos(z) - z * sqrt(omz2);
	      if(contrib >= 0.0) {
		total += contrib;
	      }
	    }
	  }
	}
      }
      values[j] = M_2_PI * total;
    }
  }
}


