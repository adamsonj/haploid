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


void
add_rec_entries (double ** sums, rtable_t ** rec_table, size_t geno)
{
  /* add the entries in the array of recombination tables rec_table;
     enter the sums into a matrix sums; this is to check that all
     entries add to 1.0 over range of the entire array (the entries in
     each table are probabilities and therefore over thesample space
     of mated pairs must add to 1) */
  for (int i = 0; i < geno; i++)
    for (int j = 0; j < geno; j++)
      {
	sums[i][j] = 0.0;
	for (int k = 0; k < geno; k++)
	  sums[i][j] += sparse_get_val (rec_table[k], i, j);
      }
  /* this will allow us to see which ones have problems */
}

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
  rec_mating (freq, &rec_test_data);

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
      /* assert (islessequal (alleles[j] - alleles_new[j], DBL_MIN)); */
    }
  /* assert (islessequal (tot, 1.0)); */

  /* now check the sums of entries in the recombination table and
     print them for inspection (later we can make this an assertion) */
  double ** sums = malloc (GENO * sizeof (double *));
  for (int i = 0; i < GENO; i++)
    sums[i] = calloc (GENO, sizeof (double));
  
  add_rec_entries (sums, rec_test_data.rec_table, GENO);
  for (int i = 0; i < GENO; i++)
    {
      for (int j = 0; j < GENO; j++)
	printf ("%9.8f ", sums[i][j]);
      printf ("\n");
    }  

  return 0;
}
