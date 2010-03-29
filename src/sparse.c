/*

  sparse.c: sparse-matrix creation, access and operations

  Copyright 2010 Joel J. Adamson 

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

#include "haploidpriv.h"

/* this file implements a compound sparse-matrix data structure
   (mainly) for representing recombination tables.

   Each such data structure is a linked list that tells accessor
   routines

   (a) the value

   (b) where it is in its "matrix" and

   (c) the location of the next non-zero element (NULL for the final
   element in a matrix); using this technique we can speed up
   multiplication by searching within this data structure, instead of
   iterating over the (non-sparse) multiplicand

   Each full recombination table will be an array of pointers to the
   first element of each matrix-representing list

*/



sparse_elt_t *
sparse_new_elt (int * indices, double value, sparse_elt_t * next)
{
  /* return a pointer to a new sparse-matrix element */
  sparse_elt_t * new_elt;
  if (( new_elt = malloc (sizeof (sparse_elt_t))) == NULL )
    /* null pointer bad bad!! exit! die! */
    error (0, ENOMEM, "Null pointer\n");

  if (indices == NULL)
    /* howdy there, need a new pointer? */
    new_elt->indices = malloc (2 * sizeof (int));
  else  
    /* set new values: */
    new_elt->indices = indices;

  /* set value: */
  new_elt->val = value;
  /* set next; may be NULL */
  new_elt->next = next;

  /* return the new element: */
  return new_elt;
}

_Bool
sparse_iselt (sparse_elt_t * list, int row, int col)
{
  /* look for elements with specified row and column in list starting
     at sparse_elt */

  /* is this function useful at all? */
  
  /* this function is recursive: what kind of performance hit do we
     get for using recursion on an O(n) procedure? */

  /* First criterion: if the row or column given is smaller than the
     row and column of the first element, then we can return right
     away */
  
  /* initialize a pointer to the first element of the list */
  int * idx = list->indices;
  /* now treat idx as an array */
  /* does the row match? */
  if (row < idx[0])
    return false;
  else if (row > idx[0])
    /* recurse */
    sparse_iselt (list->next, row, col);
  else if (col < idx[1])
    /* otherwise the row does match and we need to check the column;
       return false if col is less than idx[1] */
    return false;
  else if (col > idx[1])
    /*   recurse if it doesn't match */
    sparse_iselt (list->next, row, col);

  /* otherwise we have a match! */
  return true;
}

double
sparse_get_val (sparse_elt_t * list, int row, int col)
{
  /* get the value at (row,col) */
  for (; list != NULL; list = list->next)
    {
      if ((list->indices[0] == row) &&
	  (list->indices[1] == col))
	return list->val;
    }
  /* this is really what it is: */
  return 0.0;
}

sparse_elt_t *
sparse_mat_mat_kron (int len, double * dense[len], sparse_elt_t * sparse)
{
  /* find Kronecker product of dense matrix DENSE and sparse matrix
     SPARSE  */

  /* the algorithm is simple; produce a new sparse matrix by iterating
     over the elements of SPARSE, copying each one and multiplying the
     values from the corresponding dense matrix entries

     The Kronecker product is commutative, so the order of arguments
     is significant only in the sense that the two C types are
     distinct  */
  sparse_elt_t * result = sparse_new_elt (NULL, 0.0, NULL);
  sparse_elt_t * endptr = result;
  sparse_elt_t * sparseptr = sparse;
  int row, col;
  double newval;
  for (; sparseptr->next != NULL;
       sparseptr = sparseptr->next)
    {
      /* copy everything interesting into the endptr */
      row =  sparseptr->indices[0];
      col =  sparseptr->indices[1];
      newval = (sparseptr->val) * dense[row][col];
      /* don't save if the result is zero */
      if (newval != 0.0)
	{
	  endptr->indices[0] = row;
	  endptr->indices[1] = col;
	  endptr->val = newval;
	  endptr->next = sparse_new_elt (NULL, 0.0, NULL);
	  endptr = endptr->next;
	}
      else continue;
    }
  /* no need to free sparseptr, since it will be null */
  free(endptr);
  return result;
}

double
sparse_mat_tot (sparse_elt_t * sparse)
{
  /* total the entries of SPARSE: this is equivalent to

     1*S*1^T

     i.e. multiplying on the right by a column of ones, and
     multiplying on the left with a row of ones; this is the operation
     needed for the recombination algorithm */
  double result = 0.0;
  sparse_elt_t * endptr = sparse;
  /* iterate along SPARSE, placing a sum in result */
  for (; endptr != NULL; endptr = endptr->next)
    result += endptr->val;

  free (endptr);
  return result;
}

void
sparse_vec_to_array (sparse_elt_t * sparse, double * arr, int len)
{
  /* transform a sparse list (matrix or array) into a one-dimension
     array (vector); this should work with either input type, but you
     will only get a vector out */
  int i;
  /* iterate over the rows of SPARSE (it is treated as a column
     vector) */
  for (i = 0; i < len; i++)
    arr[i] = sparse_get_val (sparse, i, 0);
}
