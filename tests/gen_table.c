/*

  gen_table.c: testing rec_gen_table
  Copyright 2009 Joel J. Adamson 

  $Id$

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
#include "../src/haploid.h"
#include "../src/sparse.h"
#define NLOCI 2
#define GENO 4
double r = 0.25;

int
main (void)
{
  sparse_elt_t ** gtable;
  gtable = rec_gen_table (NLOCI, GENO, &r);
  /* gtable is now an array of length GENO */

  int i, j, k;
  for (k = 0; k < GENO; k++)
    {
      fprintf (stdout, "Offspring %i:\n", k);
      /* print rows: */
      for (i = 0; i < GENO; i++)
	{
	  for (j = 0; j < GENO; j++)
	    fprintf (stdout, "%9.8f  ",
		     sparse_get_val (gtable[k], i, j));
	  fprintf (stdout, "\n");
	}
    }
  
  return 0;
}

