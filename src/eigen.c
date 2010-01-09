/*

  eigen.c: calculating jacobians and eigenvalues for haploid two-sex models
  Copyright 2009 Joel J. Adamson 
  $Id: eigen.c 497 2009-06-20 00:23:00Z joel $

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

/*

  solving for the eigenvalues may require the GNU Scientific Library

*/

#include "func.h"
#include "haploid.h"
static int build_jac (Xstate_t * X, double *matrix, int rows, int cols);
static inline double dDdxi (Xstate_t * X, int i);

int
get_eigen (Xstate_t * X)
{
  /* 

     in the future, we can arrange to store the data in variable
     length arrays, that the program can determine from the number of
     input columns (number of haplotypes), but for right now, this
     matrix is 4x4

     all these 4s can be replaced later

   */
  double *matrix = malloc (16 * sizeof (double));
  if (matrix != NULL)
    {
      build_jac (X, matrix, 4, 4);
      /* 

         okay, we've built the jacobi matrix, now get its eigenvalues
         and save them for printing

       */
      gsl_matrix_view v = gsl_matrix_view_array (matrix, 4, 4);
      gsl_vector_complex *eigval = gsl_vector_complex_alloc (4);
      gsl_matrix_complex *eigvec = gsl_matrix_complex_alloc (4, 4);

      gsl_eigen_nonsymmv_workspace *wspace = gsl_eigen_nonsymmv_alloc (4);
      /* the actual calculation: */
      gsl_eigen_nonsymmv (&v.matrix, eigval, eigvec, wspace);

      gsl_eigen_nonsymmv_free (wspace);
      gsl_eigen_nonsymmv_sort (eigval, eigvec, GSL_EIGEN_SORT_ABS_DESC);

      /*

         Xmat_print is actually a generic array printer: Yippee!

         as long as these things are arrays, I can print them using Xmat_print

         Hold it, no Xmat_print assumes two-dimensional arrays


         I may want to dress it up a little, for example saying
         "Eigenvalue = ..." or "\lambda_1 = ..." or "\lambda_%s = "

       */
      Xvec_print (eigval->data, "\\lambda =");
      Xmat_print (eigvec->data, "Eigenvectors:");
      Xmat_print (matrix, "Jacobi matrix: ");

      free (matrix);
      gsl_vector_complex_free (eigval);
      gsl_matrix_complex_free (eigvec);
    }
  else
    {
      perror ("Null pointer");
      return ENOMEM;
    }

  return 0;
}

/*

  delegates calculation of jacobian and eigenvalues of jacobians given
  an Xstate X; returns 0 or error code; modifies X

  This function should only delegate and return, not calculate (more
  complex)

  this function needs to return or modify what will become the printed
  representation of the jacobian and the eigenvalues and eigenvectors

*/
static int
build_jac (Xstate_t * X, double *matrix, int rows, int cols)
{
  /*

     calculate the Jacobi matrix

     iterate over the haplotypes and populate the matrix with the
     derivative values

     see macros for derivative calculations

     dDdxi returns a pointer to an array containing values for the
     derivative of LD with respect to x_i; iterate over this pointer
     while calculating the derivatives (will be constant for columns)

   */

  int eta = 1;
  int i, j;			/* counters for rows and columns */
  for (i = 0; i < rows; i++)
    {
      if (i == 0 || i == 3)
	eta = 1;
      else
	eta = -1;

      for (j = 0; j < cols; j++)
	{
	  if (i == j)
	    {
	      *matrix = DFIDXI (X->p[i], X->w[i], X->wbar,
				eta, X->r, dDdxi (X, i));
	    }
	  else
	    {
	      *matrix = DFIDXJ (X->p[i], X->w[i], X->w[j],
				X->wbar, eta, X->r, dDdxi (X, j));
	    }
	  matrix++;
	}
    }
  return 0;
}

/*

  translate a regular array matrix to a GSL matrix

  use the gsl matrix for solving for eigen values and eigenvectors

  use the regular matrix for printing purposes

  looks like this is easily accomplished with gsl_matrix_view_array,
  so this function may be unnecessary

  this function should go ahead an return the eigenvalues and eigenvectors

*/
static inline double
dDdxi (Xstate_t * X, int i)
{
  /*

     calculate the derivative of LD with respect to x_i

   */

  double xj = 0.0;

  switch (i)
    {
      /* use cases to define the components, then add them at the end */
    case 3:
    case 0:
      {
	if (i == 3)
	  xj = X->p[0];
	else
	  xj = X->p[3];
      }
    case 2:
    case 1:
      {
	if (i == 2)
	  xj = X->p[1];
	else
	  xj = X->p[2];
      }
    }

  return xj;
}
