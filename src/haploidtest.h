#ifndef HAPLOIDTEST_H
#define HAPLOIDTEST_H

/*

  haploid.h: Public header for haploid

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

/* includes */
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <stdbool.h>
#include <limits.h>
#include <errno.h>
#include <error.h>

/* sparse.c */

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
  int geno;			/* number of genotypes */
  int nloci;			/* number of loci */
  rtable_t ** rec_table;	/* recombination table */
  double ** mtable;		/* mating table (matrix) */
};


void
sparse_vec_to_array (sparse_elt_t * sparse, double * arr, int len);

sparse_elt_t *
sparse_new_elt (int * indices, double value, sparse_elt_t * next);

double
sparse_get_val (sparse_elt_t * list, int row, int col);

sparse_elt_t *
sparse_mat_mat_kron (int len, double * dense[len], sparse_elt_t * sparse);

double
sparse_mat_tot (sparse_elt_t * sparse);
/* spec_funcs.c */
int
sim_stop_ck (double * p1, double * p2, int len, long double tol);

double
gen_mean (double * props, double * vals, int geno);

/* rec.c */
double *
rec_mating (haploid_data_t * data);

rtable_t **
rec_gen_table (int nloci, int geno, double * r);

double
rec_total (unsigned int nloci, unsigned int diff,
	   double * r, _Bool recomb_p);

/* geno_func.c */
void
allele_to_genotype (double * allele_freqs, double * geno_freqs,
		    int nloci, int geno);

void
genotype_to_allele (double * allele_freqs, double * geno_freqs,
		    int nloci, int geno);

/* mating.c */
double **
rmtable (int geno, double * freq);

/* bits.c: useful functions for integers */

_Bool
bits_isset (int x, unsigned int pos);

unsigned int
bits_extract (unsigned int start, unsigned int end, unsigned int x);

unsigned int
bits_popcount (int x);

unsigned int
bits_ffs (unsigned int x);



#endif	/*  HAPLOIDTEST_H */
