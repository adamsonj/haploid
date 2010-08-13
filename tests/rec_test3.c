/*

  rec_test.c: testing recombination algorithms
  Copyright 2009 Joel J. Adamson 

  $Id: rec_test.c 1046 2010-05-07 19:32:41Z trashbird1240 $

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
#include <float.h>
#include <assert.h>
#include <stdio.h>
#include "../src/haploid.h"
#include "../src/sparse.h"
#include <time.h>

#define NLOCI 3
#define GENO 8
#define LEN 2
#ifndef R
#define R 0.5F
#endif
double r[LEN] = { R, R};

rtable_t ** rec_table;

/* declarations: */
void
rec_test_total (void);

void
rec_test_prtable (haploid_data_t * data);

int
main (void)
{
  rec_table = rec_gen_table (NLOCI, GENO, r);
  double alleles[3] = { 0.5, 0.5, 0.5};
  
  double freq[GENO];
  allele_to_genotype (alleles, freq, NLOCI, GENO);
  
  haploid_data_t rec_test_data =
    { GENO, NLOCI, rec_gen_table (NLOCI, GENO, r), rmtable (GENO, freq)};

  /* time the next calculation */
  rec_test_prtable (&rec_test_data);

  double tot = 0.0F;
  /* time the next calculation */
  time_t time1 = time (NULL);
  rec_mating (freq, &rec_test_data);
  time_t time2 = time (NULL);
  /* print freq */
  for (int j = 0; j < GENO; j++)
    {
      printf ("x[%1x] = %9.8f\n", j, freq[j]);
      tot += freq[j];
    }
  /* print alleles */
  double alleles_new[NLOCI];
  genotype_to_allele (alleles_new, freq, NLOCI, GENO);
  for (int j = 0; j < NLOCI; j++)
    {
      printf ("p[%1x] = %9.8f\n", j, alleles_new[j]);
      assert (islessequal (alleles[j] - alleles_new[j], DBL_MIN));
    }
  assert (islessequal (tot, 1.0));

  printf ("Call to rec_mating () took %9.8f sec\n", difftime(time2, time1));
  return 0;
}

