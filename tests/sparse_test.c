/*

  sparse_test.c: testing sparse_matrix operations
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

#include "../src/haploidtest.h"


#define LEN 0xf

int
main (void)
{
  /* matrix to multiply times identity matrix */
  double ** matelts = calloc (LEN, sizeof (double *));
  if (matelts == NULL) error (0, ENOMEM, "Null pointer\n");
  
  /* identity matrix: */
  int matindx[LEN][2];
  /* array of matrix elements */
  sparse_elt_t * mat[LEN];
  int i, j;

  for (i = 0; i < LEN; i++)
    {
      for (j = 0; j < 2; j++)
	matindx[i][j] = i;

      mat[i] = sparse_new_elt (matindx[i], 1.0, NULL);
      /* set the ->next of each element, except the last one, which
	 should stay NULL */
      if (i > 0)
	mat[i-1]->next = mat[i];

      /* now fill in the matrix for multiplication */
      matelts[i] = calloc (LEN, sizeof (double));
      if (matelts[i] == NULL)
	error (0, ENOMEM, "Null pointer\n");

      for (j = 0; j < LEN; j++)
	matelts[i][j] = i;
    }

  sparse_elt_t * matresult = sparse_mat_mat_kron (LEN, matelts, *mat);
  for (sparse_elt_t * matptr = matresult;
       matptr->next != NULL;
       matptr = matptr->next)
    printf ("M[%2i,%2i] = %9.8f\n",
	    /* row index */
	    matptr->indices[0],
	    /* column index */
	    matptr->indices[1],
	    /* value */
	    matptr->val);

 
  double res;
  res = sparse_mat_tot (matresult);

  printf ("Total of entries in M = %9.8f \n", res);

  return 0;
}
	 
