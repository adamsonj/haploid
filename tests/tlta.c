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
#define R 0.0374		/* R */

/* includes */
#include "../src/haploid.h"
#include "../src/bithacks.h"

int
next_gen (double * freqs, double * W);

/* recombination table */
double rect[GENO][GENO][GENO];


int
main (void)
{

  double old[NLOCI];
  double freq[GENO];
  double allele[NLOCI];

  int i, j, stop_p, break_p, n;

  /* some initializations */

  /*    additive fitness for each genotype:
	W(ab) = 1
	W(Ab) = 1 + s_1
	W(aB) = 1 + s_2
	W(AB) = 1 + (s_1 + s_2)  */ 

  double W[GENO] = { 1.0, 1.3, 1.5, 1.8};

  /* initialize recombination table: */
  double rprob = R;
  set_rec_table (NLOCI, GENO, rect, &rprob);

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
	      if (B_IS_SET (i, j))
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
	  
      /* while sim_stop_ck returns 1 and we are at less than a million
	 generations, keep going, baby */
      while (stop_p && n < GENS)
	{

	  /* first print the allele frequencies */
	  for (j = 0; j < NLOCI; j++)
	    printf ("%16.15f ", allele[j]);
	  
	  /* print linkage disequilibrium; valid only for TLTA */
	  printf ("%16.15f ", freq[0] * freq[3] - freq[1] * freq[2]);

	  /* flush the output */
	  printf ("\n");
	  
	  /* save the old values for comparison */
	  for (j = 0; j < NLOCI; j++)
	    old[j] = allele[j];
	  
	  /* produce the next generation */
	  /* be sure to check the return value of next_gen :) */
	  break_p = next_gen (freq, W);
	  if (break_p) break;
	  /* compare old and new */
	  else
	    {
	      /* generate new allele frequencies: */
	      genotype_to_allele (allele, freq, NLOCI, GENO);
	      stop_p = sim_stop_ck (allele, old, NLOCI, 1e-8);
	      n++;
	    }
	}
    }	
  printf ("You win!\n");
  return 0;
}

int
next_gen (double * freqs, double * W)
{
  /* selection step: Multiply each genotype frequency by its
     respective fitness, then divide by the mean fitness  */

  /* new frequencies after selection step */
  double new0[GENO];
  double wbar = gen_mean (freqs, W, GENO);
  int i;
  for (i = 0; i < GENO; i++)
    new0[i] = freqs[i] * W[i] / wbar;

  /* mating step: generate mating table, then send to recombination */

  /* new frequencies after mating step */
  /* populate the arrays with zeros */
  double new[GENO] = {[0] = 0};
  /* use frequencies in freqs to generate mating table: */
  double F[GENO][GENO];
  
  /* random mating: rmtable initializes F and recombination performs
     sex*/
  rmtable (GENO, new0, F);
  recombination (GENO, rect, F, new);
      
  /* update freqs with new frequencies */
  for (i = 0; i < GENO; i++)
    freqs[i] = new[i];
  return 0;
}
