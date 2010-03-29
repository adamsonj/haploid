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

double *
selection (double * freqs, double * W);

void
rec_test_prtable (rtable_t ** rtable);


int
main (void)
{
  

  double old[NLOCI];
  double * freq = calloc (GENO, sizeof (double));
  double allele[NLOCI];

  int i, j, stop_p, n;

  /* some initializations */

  /*    additive fitness for each genotype:
	W(ab) = 1
	W(Ab) = 1 + s_1
	W(aB) = 1 + s_2
	W(AB) = 1 + (s_1 + s_2)  */ 

  double W[GENO] = { M_1_PI, M_PI_4, M_PI_2, M_PI};

  /* initialize recombination table: */
  double rprob = R;
 

  /* initialize random number generator; this is a literal since we
     want to replicate values on every trial */
  srand48 (0);
  for (i = 0; i < TRIALS; i++)
    {
      if (i < GENO)
	{
	  /* for the first GENO cases we want to test the behavior of
	     the model under the degenerate equilibria: if we have
	     only one genotype in the population

	     We need to set the allele frequencies to the
	     corresponding bits of i; e.g. for i = 0, both allele
	     frequencies should be 0; for i = 1, the first allele
	     should have frequency 0, the second should have frequency
	     1 */
	  for (j = 0; j < NLOCI; j++)
	    {
	      if (bits_isset (i, j))
		allele[j] = 1;
	      else
		allele[j] = 0;
	    }
	}
      else
	{
	  /* for the rest of the cases, assign random numbers */
	  for (j = 0; j < NLOCI; j++)
	    allele[j] = drand48 ();
	}	  
      
      /* print a header for each case */
      printf ("\n\n");

      /* preliminaries: */
      stop_p = 1;      
      n = 0;
      /* install the initial genotype frequencies to freq from allele */
      allele_to_genotype (allele, freq, NLOCI, GENO);
      haploid_data_t tlta_data =
	{GENO, NLOCI,
	 rec_gen_table(NLOCI, GENO, &rprob)};
      /* while sim_stop_ck returns 1 and we are at less than a million
	 generations, keep going, baby */
      while (stop_p && n < GENS)
	{

	  /* first print the allele frequencies */
	  for (j = 0; j < NLOCI; j++)
	    printf ("%9.8f ", allele[j]);
	  
	  /* print linkage disequilibrium; valid only for TLTA */
	  printf ("%9.8f ", freq[0] * freq[3] - freq[1] * freq[2]);

	  /* flush the output */
	  printf ("\n");
	  
	  /* save the old values for comparison */
	  for (j = 0; j < NLOCI; j++)
	    old[j] = allele[j];
	  
	  /* produce the next generation */
	  tlta_data.mtable = rmtable (GENO,  selection (freq, W));
	  freq = rec_mating (&tlta_data);
	  /* compare old and new */
	  /* generate new allele frequencies: */
	  genotype_to_allele (allele, freq, NLOCI, GENO);
	      
	  stop_p = sim_stop_ck (allele, old, NLOCI, 1e-8);
	  n++;

	}
    }	
  printf ("You win!\n");
  return 0;
}

double *
selection (double * freqs, double * W)
{
  /* selection step: Multiply each genotype frequency by its
     respective fitness, then divide by the mean fitness  */

  /* new frequencies after selection step */
  double * new = calloc (GENO, sizeof (double));
  double wbar = gen_mean (freqs, W, GENO);
  int i;
  /* selection: */
  for (i = 0; i < GENO; i++)
    new[i] = freqs[i] * W[i] / wbar;

  /* update freqs with new frequencies */
  return new;
}

void
rec_test_prtable (rtable_t ** rtable)
{
  /* traverse the table, printing zeros where there are no entries */
  sparse_elt_t * tptr;
  int geno;			/* offspring genotype */
  double val;
  int i, j;			/* row, column indices */
  for (geno = 0; geno < GENO; geno++)
    /* iterate over the rows of the table, printing the whole thing as
       a matrix (with zero entries) */
    {
      tptr = rtable[geno];

      /* print a header announcing the genotype: */
      printf ("Offspring %x:\n", geno);
      for (i = 0; i < GENO; i++)
	{
	  /* print the parent genotypes across the top */
	  for (j = -1; j < GENO; j++)
	    if (j < 0)
	      printf ("%2s", " ");
	    else
	      printf ("%10x ", j);
	  printf ("\n");
	  /* print the paternal genotype for the row */
	  printf ("%x:", i);
	  /* now print probabilities */
	  for (j = 0; j < GENO; j++)
	    {
	      val = sparse_get_val (tptr, i, j);
	      printf ("%9.8f ", val);
	    }
	  printf ("\n");
	}
      printf ("\n");
    }	  
}
