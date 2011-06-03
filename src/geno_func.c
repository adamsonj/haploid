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
#include <stdlib.h>
#include "haploid.h"

void
allele_to_genotype (double * allele_freqs, double * geno_freqs,
		    size_t nloci, size_t geno)
{
  /* take an array of allele frequencies and generate an array of
     genotype frequencies; store the result in geno_freqs */


  for (int i = 0; i < geno; i++)
    {
      double p = 1.0;
      
      for (int j = 0; j < nloci; j++)
	{
	  if (bits_isset (i, j))
	    p *= allele_freqs[j];
	  else
	    p *= fdim(1.0, allele_freqs[j]);
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
  for (int i = 0; i < nloci; i++)
    {
      double p = 0.0;
      
      for (int j = 0; j < geno; j++)
	{
	  if (bits_isset (j, i))
	    p += geno_freqs[j];
	}

      allele_freqs[i] = p;
    }
}

double
ld_from_geno (double * genofreqs, size_t geno)
{
  /* takes a single argument genofreqs and calculates the linkage
     disequilibrium by finding the difference between the first
     element and the product of the allele frequencies making up that
     least-significant genotype */
  size_t nloci = (uint) log2 (geno);
  double alleles[nloci];
  /* assume diallelic loci and calculate the genotype frequencies */
  genotype_to_allele (alleles, genofreqs, nloci, geno);
  double minuend = *genofreqs;
  double subtrahend = 1.0;
  for (int i = 0; i < nloci; i++)
    subtrahend *= fdim(1.0, alleles[i]);
  return fdim(minuend, subtrahend);
}

double
ld_sub_geno (double * genofreqs, uint loci, size_t ngeno)
{
  /* GENOFREQS is a vector of genotype frequencies; LOCI is the
     genotype (a mask) that has both alleles in the "on" position;
     NGENO is the number of elements of GENOFREQS */
  /* get "sub-genotype" frequency */
  double subgeno_freq = 0.0F;
  uint nloci = (uint) log2(ngeno);
  double alleles[nloci];
  genotype_to_allele (alleles, genofreqs, nloci, ngeno);
  for (int i = 0; i < ngeno; i++)
    if ((i & loci) == loci) subgeno_freq += genofreqs[i];
  /* now find the allele frequencies of each "on" allele at the
     requested loci */
  /* which allele frequencies do we want from ALLELES? */
  /* find all the set bits of LOCI */
    double subtrahend = 1.0F;
  for (int i = 0; i < nloci; i++)
    if (bits_isset (loci, i))
      subtrahend *= alleles[i];
  return subgeno_freq - subtrahend;  
}
