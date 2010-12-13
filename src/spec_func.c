/*

  spec_func.c: special mathematical functions for use in genetic
  simulations
  
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
#include <math.h>
double
gen_mean (double * props, double * vals, int len)
{
  /* calculate a generalized mean given an array of values and
     probabilities */     
  double mean = 0.0;
  for (int i = 0; i < len; i++)
    mean += props[i] * vals[i];
  /* return the mean */
  return mean;
}

static inline double
euclid_dist (double * array1, double * array2, int len)
{
  /* a simple Euclidean distance function */
  double diffarray[len];
  double eudiff = 0.0;
  for (int i = 0; i < len; i++)
    {
      diffarray[i] = array1[i] - array2[i];
      eudiff += diffarray[i] * diffarray[i];
    }
  return sqrt (eudiff);
}

int
sim_stop_ck (double * p1, double * p2, int len, long double tol)
{
  /*

    check the Euclidean distance between two arrays and if its value
    is smaller than tol, stop the simulation 

    Future version reserve the right to use different functions to
    assess doneness

  */

  if (euclid_dist (p1, p2, len) < tol)
    return 0;
  else
    /* return 1 to signal that distance is still large */
    return 1;
}

