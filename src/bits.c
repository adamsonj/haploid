/*

  bits.c: bit-wise abstractions
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

_Bool
bits_isset (int x, unsigned int pos)
{
  return (x & (1 << pos))?true:false;
}

unsigned int
bits_extract (unsigned int start, unsigned int end, unsigned int x)
{
  /* extract an unsigned integer between positions start (inclusive)
     to end (exclusive) from an unsigned integer X */
  unsigned int len = end - start;
  x >>= start; /* okay to modify local x */
  /* shift UINT_MAX (length of uint_max) - len positions to the
     right: */
  unsigned int maxshift = UINT_MAX >> (__WORDSIZE - len);
  return x & maxshift; /* exclude all the values above len
			  positions */
}

unsigned int
bits_popcount (int x)
{
#ifdef __GNUC__
  return __builtin_popcount (x);
#else  /* __GNUC__ */

   /* Based on HS Warren. 2003.  Hacker's Delight.  Addison-Wesley, */
   /* Reading, MA, pg 73 */

  return (((x * 0x0002000400080010ULL)				
	   & 0x1111111111111111ULL )
	  * 0x1111111111111111ULL )
    >> 60;
#endif	/* __GNUC__ */
}
  
