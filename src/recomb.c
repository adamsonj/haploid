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
  /*  calculates the probability of a given recombination pattern. */
  /* written by Mark Kirkpatrick */
  double prob = 1.0;
  int i;			/* counter */
  _Bool bset0, bset1;
  for (i = 0; i < nloci - 1; i++)
    {
      bset0 = bits_isset (mask, i);
      bset1 = bits_isset (mask, i+1);
      prob *= (bset0 == bset1)?(1.0 - r[i]):r[i];
    }
  return prob;
}

static void
zygote_genotypes (unsigned int mom,
		  unsigned int dad,
		  unsigned int mask,
		  unsigned int zyg[2],
		  int nloci)
{
  /* forms the 2 zygote genotypes resulting from a Mom, a Dad &
     particular recombination mask */
  /* written by Mark Kirkpatrick */
  unsigned int zyg1 = 0;
  unsigned int zyg2 = 0;
  int k;			/* counter */
  int powk;			/* mask for a position */
  for (k = 0; k < nloci; k++)
    {

      powk = 1 << k;
      /* does the mask call for this bit to come from the father? */
      if (bits_isset (mask, k)) 
	{
	  /* get the 2^kth bit from the father */
	  zyg1 |= (dad & powk);
	  /* do the same for the other zygote */
	  zyg2 |= (mom & powk);
	}
      else
	{
	  zyg1 |= (mom & powk);
	  zyg2 |= (dad & powk);
	}
    }
  /* create the zygotes  */
  zyg[0] = zyg1;
  zyg[1] = zyg2;
}

int
set_rec_table (int nloci, int geno,
	       double rec_table[geno][geno][geno], double * r)
{
  /* sets up table of probabilities for offspring from all possible
     matings:

     rec_table [mom's genotype][dad's genotype][zygote's genotype] */ 
  unsigned int zyg[2] = { 0, 0 };
  /* all masks above mask_lim have bitwise complements below
     1 << (nloci-1)

     Therefore we _do_ get all possible combinations of loci by
     setting the upper limit at this value.  Also remember that we do
     not use mask_lim as a mask (see for...mask below) */
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
	  /* initialize rec_table */
	  for (k = 0; k < geno; k++)
	    rec_table[i][j][k] = 0;

	  /* iterate over recombination possibilities */
	  for (mask = 0; mask < mask_lim; mask++)
	    {
	      /* the probability of producing either of the two
		 zygotes from the particular recombination event
		 specified by mask  */
	      zygote_prob = (0.5) * rec_prob (mask, nloci, r);
		    
	      /* Generate zygote genotypes (indices in the
		 rec_table) */
	      zygote_genotypes (i, j, mask, zyg, nloci);
	      /* modify the recombination table: */
	      for (k = 0; k < 2; k++)
		  rec_table[i][j][zyg[k]] += zygote_prob;
	    }
	}
    }
  return 0;
}



int
recombination (int geno, double rec_table[geno][geno][geno],
	       double mating_table[geno][geno],
	       double xt[geno])
{
  int i, j, k;
  for (k = 0; k < geno; k++)
    for (i = 0; i < geno; i++)
      for (j = 0; j < geno; j++)
	/* for the k-th entry in xt, we need to collect the k-th entry
	   of every pointer in rec_table, multiplied times the i,j-th
	   entry in mating_table */
	xt[k] += mating_table[i][j] * rec_table[i][j][k];
  return 0;
}
