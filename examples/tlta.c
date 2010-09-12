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
#ifndef R
#define R 0.5F
#endif
#define MAXBUF 256

/* includes */
#include <stdio.h>
#include "../src/haploid.h"
#include "../src/sparse.h"
#include <string.h>
#include <time.h>

char prec[] = "%9.8f ";

void
selection (double * freqs, double * W);

void
rec_test_prtable (haploid_data_t * data);


int
main (void)
{
  /* some initializations */

  /* additive fitness for each genotype */
  /* both alleles should fix */
  double W[GENO] = { M_1_PI, M_PI_4, M_PI_2, M_PI};
  double goal[NLOCI] = { 1.0F, 1.0F };

  /* initialize recombination table: */
  double rprob = 0.25;
  rtable_t ** rtable =  rec_gen_table(NLOCI, GENO, &rprob);
 
  for (int i = 0; i < TRIALS; i++)
    {
      size_t snck;
      size_t remain = MAXBUF;
      char outstr[MAXBUF];
      char * dest = outstr;
      
      double allele[NLOCI];
      haploid_data_t tlta_data = {GENO, NLOCI, rtable};
      srand48 (time (0));
      if (i < GENO)
	for (int j = 0; j < NLOCI; j++)
	  allele[j] = (double)bits_isset (i, j);
      else
	for (int j = 0; j < NLOCI; j++)
	  allele[j] = drand48 ();

      double freq[GENO];
      
      /* install the initial genotype frequencies to freq from allele */
      allele_to_genotype (allele, freq, NLOCI, GENO);

      
	/* first print the allele frequencies */
	snck = snprintf (outstr, remain, "Trial %i\n", i);
	if (snck > remain)
	  /* abort */
	  error (ENOMEM, ENOMEM,
		 "Failed write to buffer by %zu bytes", remain - snck);
	dest += snck;
	remain -= snck;
	for (int j = 0; j < NLOCI; j++,
	       dest += snck,
	       remain -= snck)
	  {
	    snck = snprintf (dest, remain, prec, allele[j]);
	    if (snck > remain)
	      /* abort */
	      error (ENOMEM, ENOMEM,
		     "Failed write to buffer by %zu bytes", remain - snck);
	  }
	snprintf (dest, remain, "\n");
	dest += 1;
	remain -= 1;
	
	/* while sim_stop_ck returns 1 and we are at less than a million
	   generations, keep going, baby */
      int n = 0;
      while (n < GENS)
	{
	  /* produce the next generation */
	  selection (freq, W);
	  tlta_data.mtable = rmtable (GENO,  freq);
	  rec_mating (freq, &tlta_data);
	  	  
	  /* generate new allele frequencies: */
	  genotype_to_allele (allele, freq, NLOCI, GENO);
#ifdef PRFREQS
	  if (i > 3)
	    {
	      for (int j = 0; j < NLOCI; j++,
		     dest += snck,
		     remain -= snck)
		{
		  snck = snprintf (dest, remain, prec, allele[j]);
		  if (snck > remain)
		    /* abort */
		    error (ENOMEM, ENOMEM,
			   "Failed write to buffer by %d bytes", remain - snck);
		}
	      snprintf (dest, remain, "\n", remain - snck);
	      dest += 1;
	      remain -= 1;

	    }
#endif  /* PRFREQS */
	  /* we expect something to fix or be lost */
	  if (!sim_stop_ck (allele, goal, NLOCI, 1e-8))
	    break;
	  else
	    n++;
	}
      
      /* print the final frequencies */
      for (int j = 0; j < NLOCI; j++,
	     dest += snck,
	     remain -= snck)
	{
	  snck = snprintf (dest, remain, prec, allele[j]);
	  if (snck > remain)
	    /* abort */
	    error (ENOMEM, ENOMEM,
		   "Failed write to buffer by %zu bytes", remain - snck);
	}
      snprintf (dest, remain, "\n");
      dest += 1;
      remain -= 1;
      {
	fprintf (stdout, outstr);
	fprintf (stdout, "\n");
      }
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
