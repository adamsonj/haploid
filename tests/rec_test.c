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
#include <signal.h>

rtable_t ** rec_table;
haploid_data_t * gdata;
/* declarations: */
#define TOL 1e-14
#define HELL SIGABRT
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
check_entries (haploid_data_t * data)
{
  size_t geno = data->geno;
  /* now check the sums of entries in the recombination table and
     print them for inspection (later we can make this an assertion) */
  double ** sums = malloc (geno * sizeof (double *));
  for (int i = 0; i < geno; i++)
    sums[i] = calloc (geno, sizeof (double));
  
  add_rec_entries (sums, data->rec_table, geno);
  for (int i = 0; i < geno; i++)
    {
      for (int j = 0; j < geno; j++)
	fprintf (stdout, "%9.8f ", sums[i][j]);
      fprintf (stdout, "\n");
    }  
  return 0;
}

void
abhandler (int sig)
{
  check_entries (gdata);
  signal (sig, SIG_DFL);
  raise (sig);
}

int
run_test (size_t nloci, double r)
{
  signal (SIGABRT, abhandler);
  size_t geno = 1 << nloci;
  double rarr[nloci];
  double alleles[nloci];
  for (int i = 0; i < nloci; i++)
    {
      rarr[i] = r;
      alleles[i] = 0.5;
    }
  rec_table = rec_gen_table (nloci, geno, rarr);
  
  double freq[geno];
  allele_to_genotype (alleles, freq, nloci, geno);
  
  haploid_data_t rec_test_data =
    { geno, nloci, rec_gen_table (nloci, geno, rarr), rmtable (geno, freq)};

#ifdef DEBUG
  fprintf (stdout, "%zu x %zu x %zu recombination table | r = %f\n", geno, geno, geno, r);
  rec_test_prtable (&rec_test_data);
#endif  /* DEBUG */

  double tot = 0.0F;
  rec_mating (freq, &rec_test_data);

  /* print freq */
#ifdef DEBUG
  fprintf (stdout, "Genotype frequencies:\n");
#endif
  for (int j = 0; j < geno; j++)
    {
#ifdef DEBUG
      fprintf (stdout, "x[%1x] = %9.8f\n", j, freq[j]);
#endif
      tot += freq[j];
    }
  /* print alleles */
  double alleles_new[nloci];
  genotype_to_allele (alleles_new, freq, nloci, geno);
#ifdef DEBUG

  fprintf (stdout, "Allele frequencies:\n");
#endif	/* DEBUG */
  for (int j = 0; j < nloci; j++)
    {
#ifdef DEBUG
      fprintf (stdout, "p[%1x] = %9.8f\n", j, alleles_new[j]);
#endif	/* DEBUG */
      if (isgreater (fabs(alleles_new[j] - alleles[j]), TOL))
	{
	  fprintf (stdout, "Change in allele frequency detected: \n"
		   "Starting frequency: p[%1x] = %9.8f\n"
		   "New frequency: p[%1x] = %9.8f\n",
		   j, alleles[j], j, alleles_new[j]);
	  gdata = &rec_test_data;
	  raise (HELL);
	}
    }
  return 0;
}

int
main (void)
{
  for (int i = 0; i < 0xf; i++)
    for (double r = 0.0F; r < 0.6; r += 0.1)
      run_test (i, r);
  return 0;
}
