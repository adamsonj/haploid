/*

  sim.c: calculations for simulations
  Copyright 2009 Joel J. Adamson 
  $Id: sim.c 497 2009-06-20 00:23:00Z joel $

  Joel J. Adamson	-- http://www.unc.edu/~adamsonj
  University of North Carolina at Chapel Hill
  CB #3280, Coker Hall
  Chapel Hill, NC 27599-3280
  <adamsonj@email.unc.edu>

  This file is part of haploid

  haploid is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  haploid is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with haploid.  If not, see <http://www.gnu.org/licenses/>.
  
*/

#include "func.h"
#include "haploid.h"

static inline double euclid_dist (double *array1, double *array2, int len);

int
sim (Xstate_t * X, int max_tries)
{
  /*

     X is an Xstate

     max_tries is the tolerable number of lines to print if the system
     fails to converge

   */
  double oldX[4] = { 0, 0, 0, 0 };
  double newX[4] = { 0, 0, 0, 0 };
  int i;			/* counter */
  int tries = 0;
  _Bool keepgoing = false;

  do
    {
      /* save the state of X */
      for (i = 0; i < 4; i++)
	oldX[i] = X->p[i];

      /* print X */
      Xst_print (X, tries);

      /* update X */
      update_X (X);

      /* analyze the new state of X */
      for (i = 0; i < 4; i++)
	newX[i] = X->p[i];
      tries++;
      /*

         keep going unless the system converges, or we have tried too
         many times

       */
      keepgoing = ((euclid_dist (oldX, newX, 4) > TOLEPS)
		   && (tries < max_tries));
    }
  while (keepgoing);

  if (tries == max_tries)
    fprintf (stdout, "\nWARNING: failure to converge.\n");

  return 0;
}
