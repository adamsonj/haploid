#ifndef HAPLOID_H
#define HAPLOID_H
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
#include "bithacks.h"

/* spec_funcs.c */
int
sim_stop_ck (double * p1, double * p2, int len, long double tol);

double
gen_mean (double * props, double * vals, int geno);

/* recomb.c */
int
recombination (int geno, double rec_table[geno][geno][geno],
	       double F[geno][geno],
	       double xt[geno]);;

int
set_rec_table (int nloci, int geno,
	       double rec_table[geno][geno][geno], double * r);
/* user must define NLOCI and GENO */

/* geno_func.c */
void
allele_to_genotype (double * allele_freqs, double * geno_freqs,
		    int nloci, int geno);

void
genotype_to_allele (double * allele_freqs, double * geno_freqs,
		    int nloci, int geno);

/* mating.c */
void
rmtable (int geno, double * freq, double table[geno][geno]);
  
#endif
