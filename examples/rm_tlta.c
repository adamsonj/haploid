/*

  rm_tlta: random mating with no selection TLTA
  
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

#include "../src/haploid.h"
#include <stdio.h>
#include <error.h>
#include <time.h>
#include <limits.h>
#include <omp.h>

#define NLOCI 2
#define GENO 4
#define TRIALS 256
#define MAXBUF 163840

static const char prec[] = "%-#9.8f ";
static const char precn[] = "%-#9.8f \n";

double r =  0.25;

char *
rm_iterate (haploid_data_t * rm_data, double * alleles, double D);

int
main (void)
{
  /* a few prelims: */
  haploid_data_t * rm_data = malloc (sizeof (haploid_data_t));
  if (rm_data == NULL)
    error (ENOMEM, ENOMEM, "Null pointer");
  rm_data->rec_table = rec_gen_table (NLOCI, GENO, &r);
  rm_data->geno = GENO;
  rm_data->nloci = NLOCI;
  srand48 (time (0));
#pragma omp parallel for
  for (int i = 0; i < TRIALS; i++)
    {
      double allele[NLOCI];
      
      for (int j = 0; j < NLOCI; j++)
	allele[j] = 0.1;
      double D;
      if (i == 0)
	D = 0;
      else
	D = drand48 () / 10.0F;

      char * output = rm_iterate (rm_data, allele, D);
#pragma omp critical(pr)
      {
	fprintf (stdout, output);
	free (output);
	fprintf (stdout, "\n");
      }
    }
  return 0;
}

char *
rm_iterate (haploid_data_t * rm_data, double * alleles, double D)
{
  size_t remain = MAXBUF * CHAR_BIT;
  char * output = malloc (remain);
  if (output == NULL)
    error (ENOMEM, ENOMEM, "Null pointer");
  char * dest = output;
  int snck;
  size_t rtot = 0;

  double genotypes[GENO];
  double old[GENO];

  allele_to_genotype (alleles, genotypes, NLOCI, GENO);
  int eta;
  for (int j = 0; j < GENO; j++)
    {
      if ((j == 0)|| (j == 3))
	eta = 1;
      else
	eta = -1;
      genotypes[j] += eta * D;
    }
  
  do
    {
      for (int j = 0; j < GENO; j++)
	old[j] = genotypes[j];
      
      rm_data->mtable = rmtable (GENO, genotypes);
      rec_mating (genotypes, rm_data);
      for (int j = 0; j < GENO; j++,
	     dest += snck,
	     remain -= snck,
	     rtot += snck)
	{
	  snck = snprintf (dest, remain, prec, genotypes[j]);
	  if (snck > remain)
	    /* abort */
	    error (ENOMEM, ENOMEM,
		   "Failed write to buffer by %d bytes", remain - snck);
	}
      snck = snprintf (dest, remain, prec,
		       ld_from_geno (NLOCI, GENO, genotypes));
      remain += snck;
      rtot += snck;
      remain -= snck;
      dest += snck;
      D *= (1 - r);
      snck = snprintf (dest, remain, precn, D);		       
      remain += snck;
      rtot += snck;
      remain -= snck;
      dest += snck;
    } while (sim_stop_ck (old, genotypes, GENO, 1e-9));
  
  return realloc (output, rtot);
}
