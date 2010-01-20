/*

  bits.c: manipulating bits and debugging integer arithmetic
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

#ifndef BUFSIZE
#define BUFSIZE sizeof (unsigned int)
#endif

char *
printbits (unsigned int n)
{

  /* returns a string with a binary representation of unsigned integer
     n */

  unsigned int i;

  char * buf = malloc (BUFSIZE * CHAR_BIT + 1);
  if (buf == NULL)
    perror ("Null pointer");
  else
    {
      for (i = 0; i < BUFSIZE + 1; i++)
	buf[i] = (B_IS_SET(n, BUFSIZE - i)) ? '1' : '0';
  
      buf[BUFSIZE + 1] = '\0';
      return buf;
    }
}
