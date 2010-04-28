/*
  
  tlta.c: Two locus two allele test simulation
  Copyright 2009 Joel J. Adamson 

  $Id$

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

/* declarations */

/* default values are for a two-locus-two-allele example */
#define NLOCI 2			/* NLOCI */
#define GENO 4			/* GENO */
#define TRIALS 10		/* TRIALS */
#define GENS 1e6		/* GENS */
#define R 0.25			/* R */

/* includes */
#include "../src/haploidtest.h"
#include <string.h>
#include <time.h>
#include <omp.h>

char prec[] = "%9.8f ";

void
selection (double * freqs, double * W);

void
rec_test_prtable (haploid_data_t * data);


int
main (void)
{
  double allele[NLOCI];

  int i, j;

  /* some initializations */

  /* additive fitness for each genotype */ 
  double W[GENO] = { M_1_PI, M_PI_4, M_PI_2, M_PI};

  /* initialize recombination table: */
  double rprob = 0.25;
  rtable_t ** rtable =  rec_gen_table(NLOCI, GENO, &rprob);
  haploid_data_t tlta_data = {GENO, NLOCI, rtable};
  srand48 (time (0));
#ifdef DEBUG
  rec_test_prtable (&tlta_data);
  exit (EXIT_SUCCESS);
#endif  /* DEBUG */
#pragma omp parallel for private (i, j, allele) firstprivate (tlta_data)
  for (i = 0; i < TRIALS; i++)
    {
      
      if (i < GENO)
	for (j = 0; j < NLOCI; j++)
	  allele[j] = (double)bits_isset (i, j);
      else
	for (j = 0; j < NLOCI; j++)
	  allele[j] = drand48 ();

      double freq[GENO];
      double old[NLOCI];
      
      /* install the initial genotype frequencies to freq from allele */
      allele_to_genotype (allele, freq, NLOCI, GENO);

#ifdef _OPENMP
      /* if we are OPENMP land, print the thread id: */
      printf ("OMP Thread %i\n", omp_get_thread_num ());
#endif	/* _OPENMP */
      
      /* first print the allele frequencies */
      printf ("Trial %i\n", i);
      for (j = 0; j < NLOCI; j++)
	printf (prec, allele[j]);
      /* flush the output */
      printf ("\n");
      /* while sim_stop_ck returns 1 and we are at less than a million
	 generations, keep going, baby */

      int n = 0;
      while (n < GENS)
	{
	  /* produce the next generation */
	  memmove (old, allele, NLOCI * sizeof (double));
	  selection (freq, W);
	  tlta_data.mtable = rmtable (GENO,  freq);
	  rec_mating (freq, &tlta_data);
	  
	  /* generate new allele frequencies: */
	  genotype_to_allele (allele, freq, NLOCI, GENO);
#ifdef PRFREQS
	  for (j = 0; j < NLOCI; j++)
	    printf (prec, allele[j]);
	  printf ("\n");
#endif  /* PRFREQS */
	  /* we expect something to fix or be lost */
	  if (!sim_stop_ck (allele, old, NLOCI, 1e-8))
	    break;
	  else
	    n++;
	}
      
      /* print the final frequencies */
      for (j = 0; j < NLOCI; j++)
	printf (prec, allele[j]);
   
      /* flush the output */
      printf ("\n");
      /* flush the toilet: */
    }
  return 0;
}

void
selection (double * freqs, double * W)
{
  /* selection step: Multiply each genotype frequency by its
     respective fitness, then divide by the mean fitness  */

  /* new frequencies after selection step */
  double wbar = gen_mean (freqs, W, GENO);
  /* selection: */
  for (int i = 0; i < GENO; i++)
    /* update freqs with new frequencies */
    freqs[i] = freqs[i] * W[i] / wbar;
}
