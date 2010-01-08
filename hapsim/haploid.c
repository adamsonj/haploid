/* Haploid Selection models

   Copyright 2009 Joel J. Adamson 
   $Id: haploid.c 482 2009-06-17 15:21:50Z joel $

   Joel J. Adamson	-- http://www.unc.edu/~adamsonj
   University of North Carolina at Chapel Hill
   CB #3280, Coker Hall
   Chapel Hill, NC 27599-3280
   <adamsonj@email.unc.edu>
   

   This is the main file containing driver routines and accessory
   subroutines

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
#include "func.h"
#include "haploid.h"

int
main (int argc, char **argv)
{
  /* print a cheerful message: */

  fprintf (stderr, "Haploid Population Genetic Numerical Analysis\n"
	   "Compiled on %s\n"
	   "Copyright 2009 Joel J. Adamson <adamsonj@email.unc.edu>\n"
	   "This is FREE SOFTWARE provided under the GNU General Public License v.3\n"
	   "This software is provided AS IS with NO WARRANTY\n"
	   "\n", __DATE__);
  io_init (argc, argv);
  return 0;
}
