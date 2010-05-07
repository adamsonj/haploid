#ifndef SPARSE_H
#define SPARSE_H

/*

  sparse.h: Public header for haploid

  Copyright 2009 Joel J. Adamson 
  
  $Id$$

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

#include "haploid.h"

sparse_elt_t *
sparse_new_elt (int * indices, double value, sparse_elt_t * next);

double
sparse_get_val (sparse_elt_t * list, int row, int col);

sparse_elt_t *
sparse_mat_mat_kron (size_t len, double * dense[len], sparse_elt_t * sparse);

double
sparse_mat_tot (sparse_elt_t * sparse);

#endif	/*  SPARSE_H */
