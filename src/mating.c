/*

  mating.c: functions to create mating tables
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

#include <stdlib.h>
#include <error.h>
#include <errno.h>
#include <assert.h>
#include <math.h>

double **
rmtable (double * freq, size_t geno)
{
  /* random mating table */
  double ** table = malloc (geno * sizeof (double *));
  if (table == NULL)
    error (0, ENOMEM, "Null pointer\n");
  double denom = 0.0F;
  for (int i = 0; i < geno; i++)
    {
      table[i] = calloc (geno, sizeof (double));
      if (table[i] == NULL)
	error (0, ENOMEM, "Null pointer\n");
      for (int j = 0; j < geno; j++)
	{
	  table[i][j] = freq[i] * freq[j];
	  denom += table[i][j];
	}
    }
  assert (isgreater (denom, 0.0));
  for (int i = 0; i < geno; i++)
    for (int j = 0; j < geno; j++)
      table[i][j] /= denom;
  return table;
}
