/*

  prtable.c: function for printing recombination tables
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
#include "../src/haploid.h"
#include "../src/sparse.h"
/* print a recombination table */

void
rec_test_prtable (haploid_data_t * data)
{
  /* traverse the table, printing zeros where there are no entries */
  rtable_t * tptr;
  rtable_t ** rtable = data->rec_table;
  size_t geno = data->geno;		/* offspring genotype */
  double val;
  int i, j, k;			/* row, column indices */
  for (k = 0; k < geno; k++)
    /* iterate over the rows of the table, printing the whole thing as
       a matrix (with zero entries) */
    {
      printf ("\n");
      tptr = rtable[k];

      /* print a header announcing the genotype: */
      printf ("Offspring %x:", k);
      /* print the parent genotypes across the top */
      for (j = -1; j < geno; j++)
	if (j < 0)
	  printf ("%2s", " ");
	else
	  printf ("%10x ", j);
      printf ("\n");
      for (i = 0; i < geno; i++)
	{
	  /* print the paternal genotype for the row */
	  printf ("%x:", i);
	  /* now print probabilities */
	  for (j = 0; j < geno; j++)
	    {
	      val = sparse_get_val (tptr, i, j);
	      printf ("%9.8f ", val);
	    }
	  printf ("\n");
	}
    }	  
}
