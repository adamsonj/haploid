/*
  trans.c: functions for numerical analysis of haploid genetics
  Copyright 2009 Joel J. Adamson 
  
  $Id: trans.c 497 2009-06-20 00:23:00Z joel $

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


  Functions for numerical analysis of haploid genetics.  This file
  contains definitions of data structures and basic functions for
  manipulating them in (a) numerical recursions and (b) numerical
  differentiations (e.g. for calculating the Jacobi matrices of the
  transformations at equilibria).
  
  The basic recurrence relation references a structure x whose members
  are the frequency (freq) of the four haplotypes, the sex-specific
  fitnesses of the four haplotypes (wm) and the average fitnesses of
  the two sexes
  
  x also contains eta = [1,-1,-1,1],  r and D
  
  So basically, the structure x contains everything needed for the
  recurrence relation
  
  The basic operation of the program will be to update x
  
*/
#include "func.h"
#include "haploid.h"

/* destructive functions */
static int update_D (Xstate_t * X);
static int update_wbar (Xstate_t * X);
static int update_p (Xstate_t * X);

static inline double D (Xstate_t * X);
static inline double wbar (Xstate_t * X, int i);
static inline int gen_mean (double *props, double *vals, int len,
			    double *result);

int
calc_results (Xstate_t * X, int cmode)
{

  /* or make this part of go_within */
  /* the results must be passed b ack to go_within for printing */
  if (cmode == ONLYSIM)
    sim (X, MAXTRIES);
  else if (cmode == ONLYEIGEN)
    get_eigen (X);
  else if (cmode == BOTH)
    {
      sim (X, MAXTRIES);
      get_eigen (X);
    }
  else
    {
      errno = ECALCMODE;
      perror ("unknown calculation mode");
    }

  return 0;
}

int
finalize_Xstate (Xstate_t * X)
{
  /*

     Check an Xstate for completeness and flag it as ready

   */

  /* local variables */
  int Disright, r_is_nan, wbar_ck;
  /* best way to  check for nans? */
  Disright = (D (X) == X->D);
  r_is_nan = isnan (X->r);


  /*

     check wbar

     should not be zero

     if it is meant to be zero, then it will not be checked again

   */
  int i;
  for (i = 0; i < 2; i++)
    {
      if (X->wbar[i] == 0)
	{
	  update_wbar (X);
	  wbar_ck = 1;
	}
      else
	wbar_ck = 1;
    }

  if (r_is_nan)
    {
      char *errfinx = malloc (50 * sizeof (char));
      sprintf (errfinx, "r is not a number: %g", X->r);
      errno = ERANGE;
      perror (errfinx);
      return ERANGE;
    }
  else if (Disright && wbar_ck)
    X->ready_flag = 1;
  else
    {
      /* update the values and confirm it ready */
      update_wbar (X);
      update_D (X);
      X->ready_flag = 1;
    }
  return 0;
}

inline double
D (Xstate_t * X)
{
  /* calculate D from the values in X */
  /*

     WARNING this function does NOT update X->D

     for that you must use the destructive function update_D

     this function may have to raise errors since it returns double

     remember to code this so you can input it as an argument to the
     gsl_function functions

   */

  double x1 = (X->p)[0];
  double x2 = (X->p)[1];
  double x3 = (X->p)[2];
  double x4 = (X->p)[3];

  double result = 0.0;
  result = (x1 * x4 - x2 * x3);
  return result;
}

inline double
wbar (Xstate_t * X, int i)
{

  /* calculate the average fitness */
  double result = 0.0;

  double w[4] = { 0, 0, 0, 0 };
  double p[4] = { 0, 0, 0, 0 };

  int j;
  for (j = 0; j < 4; j++)
    {
      w[j] = X->w[i][j];
      p[j] = X->p[j];
    }
  /*

     X->p is not a pointer, it's an array

     pass the pointer

   */

  gen_mean (p, w, 4, &result);

  return result;
}

int
update_p (Xstate_t * X)
{
  /* the basic recurrence relation */

  /* declare a few arrays: */
  double x[3] = { 0.0, 0.0, 0.0 };
  /* only need the first three entries */


  /* some scalars: */
  double w = 1;
  double wbarm = X->wbar;
  double r = X->r;
  int eta = 1;

  /* so here we go: */
  double Diseq = X->D;
  int i;
  double sum = 0, newp = 0;
  double sex_sum = 0;
  for (i = 0; i < 3; i++)
    {
      /*

         make sure not to update wbarm and wbarf (which are functions
         of x)

       */
      w = (X->w)[0][i];
      x[i] = (X->p)[i];
      if (i == 0 || i == 3)
	eta = 1;
      else
	eta = -1;

      frac = w / wbarm;
      newp = x[i] * frac - eta * r * Diseq;
      (X->p)[i] = newp;
      sum += newp;
    }

  (X->p)[3] = 1 - sum;
  return 0;
}

int
update_D (Xstate_t * X)
{
  /* update the value of D, the coefficient of linkage
     disequilibrium */

  /* call D to get the new value, and assign it to X->D */
  double newD = D (X);
  X->D = newD;
  return 0;
}

int
update_wbar (Xstate_t * X)
{
  /* update wbar */
  int i;
  for (i = 0; i < 2; i++)
    (X->wbar)[i] = wbar (X, i);
  return 0;
}

void
update_X (Xstate_t * X)
{
  /*

     this comes first, to initialize X for the haplotype
     transformations

     update_wbar requires no more information than the fitnesses of the
     two sexes

   */

  update_wbar (X);
  /*

     also to initialize X

   */
  update_D (X);
  /* finally, update with the new haplotype frequencies */
  update_p (X);
}

static inline int
gen_mean (double *props, double *vals, int len, double *result)
{
  /* calculate a generalized mean given an array of values and
     probabilities */
  int i;
  /*

     calculate the length of props and vals; raise an error if

     a. length(props) != length(vals)
     b. result is an array greater than the size of a single double

   */
  double mean = 0.0;
  for (i = 0; i < len; i++)
    mean += (*props++) * (*vals++);
  /* return the mean */
  (*result) = mean;
  /* reserve for errors */
  return 0;
}
