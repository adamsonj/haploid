/*
  sim_stop.c: test simulation stop feature in sim_stop_ck ()
  Copyright 2009 Joel J. Adamson 

  $Id$

  Joel J. Adamson	-- http://www.unc.edu/~adamsonj
  University of North Carolina at Chapel Hill
  CB #3280, Coker Hall
  Chapel Hill, NC 27599-3280
  <adamsonj@email.unc.edu>

  This file is part of haploid

  haploid is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation, either version 3 of the License, or (at your
  option) any later version.

  haploid is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
  for more details.

  You should have received a copy of the GNU General Public License
  along with haploid.  If not, see <http://www.gnu.org/licenses/>.


*/
#include <stdio.h>
#include <assert.h>
#include "../src/haploid.h"
#include "../src/sparse.h"

#define GENO 4
#define LEN 2

int
main (void)
{
  /* produce a small, silly simulation that will stop after a few
     rounds  */
  
  /* define an array */
  double simarr[LEN] = {1.0, 0.5};
  double oldarr[LEN] = {0.0, 0.0};  
  int n = 2;

  while (sim_stop_ck (simarr, oldarr, LEN, 1e-16) && n < 1e6)
    {

      /* while you reduce the size of both elements to zero, check
	 them successively using sim_stop_ck; use large changes at
	 first, then reduce them to zero */
      for (int i = 0; i < LEN; i++)
	{
	  oldarr[i] = simarr[i];
	  simarr[i] = simarr[i] / (n * n);
	}
      n++;
    }
  /* now the "manual" check: basically mimic the Euclidean distance
     function and check it against the results of sim_stop_ck */
  double diffs;
  for (int i = 0; i < LEN; i++)
    {
      double diff1 = simarr[i] - oldarr[i];
      diffs += diff1 * diff1;
    }

  assert ((n < 1e6) && (sqrt (diffs) < 1e-16));
  return 0;
}
