#include <R.h>
#include <R_ext/Utils.h>
#include "chunkloop.h"
#include "looptest.h"

/*

  Efiksel.c

  $Revision: 1.6 $     $Date: 2022/10/22 10:09:51 $

  C implementation of 'eval' for Fiksel interaction (non-hardcore part)

  Assumes point patterns are sorted in increasing order of x coordinate

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

*/

double sqrt(double x);
double exp(double x);

void Efiksel(
  /* inputs */
  int *nnsource, 
  double *xsource,
  double *ysource,
  int *nntarget,
  double *xtarget,
  double *ytarget,
  double *rrmax,
  double *kkappa,
  /* output */
  double *values
) {
  int nsource, ntarget, maxchunk, j, i, ileft;
  double xsourcej, ysourcej, xleft, dx, dy, dx2, d2;
  double rmax, r2max, r2maxpluseps, kappa, total;

  nsource = *nnsource;
  ntarget = *nntarget;
  rmax = *rrmax;
  kappa = *kkappa;

  if(nsource == 0 || ntarget == 0) 
    return;

  r2max = rmax * rmax;
  r2maxpluseps = r2max + EPSILON(r2max);

  ileft = 0;

  OUTERCHUNKLOOP(j, nsource, maxchunk, 16384) {
    R_CheckUserInterrupt();
    INNERCHUNKLOOP(j, nsource, maxchunk, 16384) {
      total = 0;
      xsourcej = xsource[j];
      ysourcej = ysource[j];
      /* 
	 adjust starting point
      */
      xleft  = xsourcej - rmax;
      while((xtarget[ileft] < xleft) && (ileft+1 < ntarget))
	++ileft;

      /* 
	 process from ileft until dx > rmax
      */
      for(i=ileft; i < ntarget; i++) {
	/* squared interpoint distance */
	dx = xtarget[i] - xsourcej;
	dx2 = dx * dx;
	if(dx2 > r2maxpluseps)
	  break;
	dy = ytarget[i] - ysourcej;
	d2 = dx2 + dy * dy;
	if(d2 <= r2max)
	  total += exp(- kappa * sqrt(d2));
      }
      values[j] = total;
    }
  }
}



