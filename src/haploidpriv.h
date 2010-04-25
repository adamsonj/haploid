#ifndef HAPLOIDPRIV_H
#define HAPLOIDPRIV_H
/*
  haploidpriv.h: functions for numerical analysis of haploid genetics
  Copyright 2009 Joel J. Adamson 
  
  $Id: func.h 497 2009-06-20 00:23:00Z joel $$

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

/* includes: */
#include <config.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <limits.h>
#include <errno.h>
#include <error.h>
#include <omp.h>


/* data types */
typedef struct sparse_elt_t sparse_elt_t;
struct sparse_elt_t {
  int * indices;		/* row, column */
  double val;			/* the value of the matrix element */
  sparse_elt_t * next;		/* pointer to next */
};

typedef sparse_elt_t rtable_t;
typedef struct haploid_data_t haploid_data_t;
struct haploid_data_t
{
  size_t geno;			/* number of genotypes */
  size_t nloci;			/* number of loci */
  rtable_t ** rec_table;	/* recombination table */
  double ** mtable;		/* mating table (matrix) */
};

/* sparse matrix operations */
void
sparse_vec_to_array (sparse_elt_t * sparse, double * arr, int len);

sparse_elt_t *
sparse_new_elt (int * indices, double value, sparse_elt_t * next);

sparse_elt_t *
sparse_mat_mat_kron (size_t len, double * dense[len], sparse_elt_t * sparse);

double
sparse_mat_tot (sparse_elt_t * sparse);
double
sparse_get_val (sparse_elt_t * list, int row, int col);
/* bitwise operations: */

_Bool
bits_isset (int x, unsigned int pos);

unsigned int
bits_ffs (unsigned int x);

unsigned int
bits_extract (unsigned int start, unsigned int end, unsigned int x);

#endif /* HAPLOIDPRIV_H */
