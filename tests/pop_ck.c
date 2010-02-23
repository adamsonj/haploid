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
  printf ("i/d |i/x |pop(i) |iso (i,4,2) |Fifth bit set?\n");
  for (i = 0; i < 0xff; i++)
    /* print the data on the integer i */
    printf ("%-3i |%02x  |%-5d  |%-3x         |%i\n", i, i,
	    bits_popcount (i), bits_extract (4, 4, i),
	    bits_isset (i, 4));
  return 0;
}
