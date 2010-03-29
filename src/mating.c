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

#include "haploidpriv.h"

double **
rmtable (int geno, double * freq)
{
  /* random mating table for haploid monoecious organisms
     (hermaphrodites) */
  int i,j;
  double ** table = malloc (geno * sizeof (double *));
  if (table == NULL)
    error (0, ENOMEM, "Null pointer\n");
  for (i = 0; i < geno; i++)
    {
      table[i] = malloc (geno * sizeof (double));
      for (j = 0; j < geno; j++)
	table[i][j] = freq[i] * freq[j];
    }
  return table;
}
