/*

  rec_test.c: testing recombination algorithms
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

#include "../src/haploidtest.h"
#include <time.h>

#define NLOCI 2
#define GENO 4
#define LEN 1
double r[LEN] = { 0.25};

rtable_t ** rec_table;

/* declarations: */
void
rec_test_total (void);

void
rec_test_prtable (haploid_data_t * data);

int
main (void)
{
  rec_test_total ();
  rec_table = rec_gen_table (NLOCI, GENO, r);
  double alleles[2] = { 0.5, 0.5};
  
  double freq[GENO];
  allele_to_genotype (alleles, freq, NLOCI, GENO);
  
  haploid_data_t rec_test_data =
    { GENO, NLOCI, rec_gen_table (NLOCI, GENO, r), rmtable (GENO, freq)};
  rec_test_prtable (&rec_test_data);

  /* time the next calculation */
  rec_test_prtable (&rec_test_data);

  /* time the next calculation */
  time_t time1, time2;
  time1 = time (NULL);
  rec_mating (&rec_test_data);
  time2 = time (NULL);
  /* print freq */
  for (int j = 0; j < GENO; j++)
    printf ("x[%1x] = %9.8f\n", j, freq[j]);

  printf ("Call to rec_mating () took %9.8f sec\n", difftime(time1, time2));
  return 0;
}

void
rec_test_total (void)
{
  printf ("Recombination probabilities: ");
  for (int i = 0; i < LEN; i++)
    printf ("%9.8f ", r[i]);
  printf ("\n");
  printf ("Probability of recombinant offspring 0x1 from parent 0xF: %9.8f\n",
	  0.5 * rec_total (NLOCI, (0x0 ^ 0x1), r, true));
  printf ("Probability of non-recombinant offspring 0x1 \n"
	  "from parents 0xF X 0x1: %9.8f\n",
	  0.5 * rec_total (NLOCI, (0xF ^ 0x1), r, false));
  printf ("Rock 'n roll dynamite, baby!.\n");
}

