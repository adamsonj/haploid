/*

  nrm.c: testing mating procedures with non-random mating
  Copyright 2009 Joel J. Adamson 

  $Id$

  Joel J. Adamson	-- http://www.unc.edu/~adamsonj
  University of North Carolina at Chapel Hill
  CB #3280, Coker Hall
  Chapel Hill, NC 27599-3280
  <adamsonj@email.unc.edu>

  This file is part of Haploid

  Haploid is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  Haploid is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with Haploid.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "../src/haploid.h"
#include <error.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>

#define GENO 4
#define TRIALS 256
#define MAXBUF 163840
static const double err = 0.1;
static const size_t nloci = 2;
static const size_t geno = GENO;
char prec[] = "%-#9.8f ";


char *
nrm_iterate (double * freqs, haploid_data_t * data);

int
main (void)
{
  double r = 0.25;
  srand48 (time (0));

#pragma omp parallel for
  for (int i = 0; i < TRIALS; i++)
    {
      haploid_data_t * nrm_data = malloc (sizeof (haploid_data_t));
      if (nrm_data == NULL)
	error (ENOMEM, ENOMEM, "Null pointer");
  
      /* populate the structure */
      nrm_data->geno = geno;
      nrm_data->nloci = nloci;
      nrm_data->rec_table = rec_gen_table (nloci, geno, &r);
  
      /* generate an assortative mating table: */
      double ** mtable = malloc (geno * sizeof (double *));
      if (mtable == NULL)
	error (ENOMEM, ENOMEM, "Null pointer");
      for (int j = 0; j < geno; j++)
	{
	  mtable[j] = malloc (geno * sizeof (double));
	  if (mtable[j] == NULL)
	    error (ENOMEM, ENOMEM, "Null pointer");
 
	}
      nrm_data->mtable = mtable;
      /* generate some frequencies */
      double freqs[geno];
      double alleles[nloci];
      for (int j = 0; j < nloci; j++)
	alleles[j] = drand48 ();
      /* transfer alleles to genotypes through one generation of random
	 mating */
      allele_to_genotype (alleles, freqs, nloci, geno);
      /* print genotype frequencies and LD */
      char * output;
      output = nrm_iterate (freqs, nrm_data);
#pragma omp critical(pr)
      fprintf (stdout, output);
      free (output);
    }
  return 0;
}

char *
nrm_iterate (double * freqs, haploid_data_t * data)
{
  /* run the simulation: return a buffer of output data */
  /* unpack the data: */
  size_t geno = data->geno;
  size_t nloci = data->nloci;
  double ** mtable = data->mtable;
  /* eventually we should have total assortative mating */
  double old[geno];

  /* how much to distort mating probability: */
  double factor;

  /* the size of the output buffer (return value of this function) */
  size_t remain = MAXBUF * CHAR_BIT;
  /* the return value (a buffer) */
  char * finstr = malloc (remain);
  if (finstr == NULL)
    error (ENOMEM, ENOMEM, "Null pointer");
  /* point to this so we can print at the right spot in the buffer */
  char * dest = finstr;
  /* keep a running total so we can realloc at the end */
  size_t snck = 0;
  size_t rtot = 0;

  /* do it! */
  do
    {

      for (int j = 0; j < geno; j++)
	{
	  /* save old frequencies */
	  old[j] = freqs[j];
	  /* adjust mating table */
	  for (int i = 0; i < geno; i++)
	    {
	      if ((bits_isset(i,0) == bits_isset(j,0))
		  && bits_isset(i,1) == bits_isset(j,1))
		factor = 1.0 - err;
	      else
		factor = err;
	      
	      mtable[i][j] = (freqs[i]*freqs[j]) * factor;
	    }
	}

      /* mating */
      rec_mating (freqs, data);

      /* print the genotype frequencies: */
      for (int j = 0; j < geno; j++,
	     dest += snck,
	     remain -= snck,
	     rtot += snck)
	{
	  snck = snprintf (dest, remain, prec, freqs[j]);
	  if (snck > remain)
	    /* abort */
	    error (ENOMEM, ENOMEM,
		   "Failed write to buffer by %d bytes", remain - snck);
	}

      /* print LD and a newline */    
      snck = snprintf (dest, remain, "%-#9.8f\n",
		       ld_from_geno (nloci, geno, freqs));
      dest += snck;
      remain -= snck;
      rtot += snck;
    } while (sim_stop_ck (old, freqs, geno, 1e-9));

  /* add a couple of newlines at the end */
  snprintf (dest, remain, "\n\n");
  dest += 2;
  remain -= 2;
  rtot += snck + 1;		/* don't forget the \0! */

  return realloc (finstr, rtot);
}
