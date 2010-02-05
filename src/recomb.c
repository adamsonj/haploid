/*
  recomb.c: recombination-related algorithms
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

/* declarations */
#include "haploidpriv.h"

static double
rec_prob (unsigned int mask, int nloci, double * r)
{
  /*  rec_prob calculates the probability of a given recombination pattern. */
  /*  The pairwise recombination rates r[0],r[1],...,r[nLoci-1] are
      declared and initialized external to this routine.  */

  /* written by Mark Kirkpatrick */

  double prob;
  int i;

  prob = 1.0;

  for (i = 0; i < nloci - 1; i++)
    prob *= (((B_IS_SET (mask, i)) == (B_IS_SET (mask, i + 1))))?
      (1.0 - r[i]):r[i];
  return prob;
}

static void
zygote_genotypes (unsigned int i,
		  unsigned int j,
		  unsigned int mask,
		  unsigned int zyg[2],
		  int nloci)
{
  
  /* zygote_genotypes */

  /* forms the 2 zygote genotypes resulting from a Mom, a Dad &
     particular recombination mask */

  /* written by Mark Kirkpatrick */

  unsigned int zyg1 = 0, zyg2 = 0;
  unsigned int bit[nloci];
  int k;
  int powk;

  /* preliminaries         */
  for (k = 0; k < nloci; k++)
    {

      powk = 1 << k;

      /* does the mask call for this bit to come from the mother? */
      bit[k] = mask & powk;

      /* Yes */
      if (bit[k] == 0)
	{
	  /*
	    If mask and powk have no 1s in common, then we get the
	    2^kth bit from the mother
	  */
	  zyg1 = zyg1 | (i & powk);
	  /*
	    then we do the same for the other zygote, getting the
	    2^kth bit from the dad
	  */
	  zyg2 = zyg2 | (j & powk);
	}
      /* No */
      else
	{
	  zyg1 = zyg1 | (j & powk);
	  zyg2 = zyg2 | (i & powk);
	}
    }
  /* creates the zygotes  */
  zyg[0] = zyg1;
  zyg[1] = zyg2;
}

int
set_rec_table (int nloci, int geno,
	       double rec_table[geno][geno][geno], double * r)
{

  /*    This 3 dimensional array sets up the recombination table of
	probabilities of the production of zygote genotypes given the
	parental genotypes.
  */

  /* rec_table [mom's genotype, dad's genotype, zygote's genotype] */ 
  unsigned int zyg[2] = { 0, 0 };
  const unsigned int mask_lim = (1 << (nloci - 1));
  double zygote_prob;

  /* counters for iteration */
  unsigned int mask;
  int i, j, k;
  /* iterate over fathers */
  for (i = 0; i < geno; i++)
    {
      /* iterate over mothers */
      for (j = 0; j < geno; j++)
	{

	  /* iterate over recombination possibilities */
	  for (mask = 0; mask < mask_lim; mask++)
	    {
	      /* probability of recombination at a particular site;
		 mask corresponds to a set of recombination locations;
		 r is the recombination map  */

	      /* the probability of producing either of the two
		 zygotes from the particular recombination event
		 specified by mask  */

	      zygote_prob = (0.5) * rec_prob (mask, nloci, r);
		    
	      /* Generate zygote genotypes
		 
		 We should get two out of using any particular mask:
		 the one from specified directly by mask, and the one
		 specified by its its bitwise complement; store the
		 result in zyg  */

	      zygote_genotypes (i, j, mask, zyg, nloci);

	      /* Find indices corresponding to the zygotes produced by
		 mask and parents, and add to their probabilities, as
		 we just found a new way to generate them  */

	      for (k = 0; k < 2; k++)
		rec_table[i][j][zyg[k]] += zygote_prob;
	    }
	}
    }
  return 0;
}



int
recombination (int geno, double rec_table[geno][geno][geno],
	       double F[geno][geno],
	       double xt[geno])
{

  /* recombination */
  /*

    takes arguments of recombination table and mating table

  */
  int i, j, k;
  for (k = 0; k < geno; k++)
    for (i = 0; i < geno; i++)
      for (j = 0; j < geno; j++)
	/* for the k-th entry in xt, we need to collect the k-th entry
	   of every pointer in rec_table, multiplied times the i,j-th
	   entry in F */
	xt[k] += F[i][j] * rec_table[i][j][k];
  return 0;
}

/*

  Erroneous output:

  00011 X 00011:0.125  0.000  0.000  1.000  

 */
