/*

  geno_func.c: genetic calculations
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

/* includes */
#include "haploidpriv.h"

void
allele_to_genotype (double * allele_freqs, double * geno_freqs,
		    size_t nloci, size_t geno)
{
  /* take an array of allele frequencies and generate an array of
     genotype frequencies; store the result in geno_freqs */
  int i, j;
  double p;

  for (i = 0; i < geno; i++)
    {
      p = 1.0;
      
      for (j = 0; j < nloci; j++)
	{
	  if (bits_isset (i, j))
	    p *= allele_freqs[j];
	  else
	    p *= (1 - allele_freqs[j]);
	}
      /* how can we ensure that we end up with frequencies that sum to
	 one? */
      /* divide by the total */
      geno_freqs[i] = p;
    }
}

void
genotype_to_allele (double * allele_freqs, double * geno_freqs,
		    size_t nloci, size_t geno)
{
  /* take an array of genotype frequencies geno_freqs and convert them
     to allele frequencies, store the result in allele_freqs */
  /* loop over alleles (loci), iterating over genotypes within: if the
     bit for that locus is set, then multiply the allele frequency for
     the locus by the genotype frequency */

  /* notice that this is the same as allele_to_genotype except for a
     few crucial differences ... ;) */
  int i, j;
  double p;
  for (i = 0; i < nloci; i++)
    {
      p = 0.0;
      
      for (j = 0; j < geno; j++)
	{
	  if (bits_isset (j, i))
	    p += geno_freqs[j];
	}

      allele_freqs[i] = p;
    }
}


