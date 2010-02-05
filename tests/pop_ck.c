/*

  pop_ck: tests bit-oriented operations
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

#include "../src/haploid.h"


int
main (void)
{
  unsigned int i;
  /* First print a header row: */
  printf ("i   |i/b   |POP(i) | ISO (i,2,2)\n");
  for (i = 0; i < 0xF; i++)
    /* print the data on the integer i */
    printf ("%-3i |%-3s | %-5d | %-3s\n", i, debug_printbits (i),
	    POP (i), debug_printbits(ISO (i, 2, 2)));
  return 0;
}
