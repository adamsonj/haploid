/*
  rec.c: recombination-related algorithms
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
#include "sparse.h"
#include <assert.h>

double *
rec_extend_r (size_t nloci, double * r)
{
  /* return a pointer to an array of values describing probabilities
     of recombination between sets I and J as described in Bürger
     (2000): I is a non-void partition of the genome including the
     first locus; J is its complement */
  /* a boundary for the iteration */
  const unsigned int geno_mask = (1 << nloci) - 1;
  double * result = malloc ((size_t) floor(geno_mask / 2) * sizeof (double));
  if (result == NULL)
    error (0, ENOMEM, "Null pointer\n");
  /* we will access the values of result by pointer arithmetic (to
     avoid calculating indices) */
  double * resultptr = result;

  /* iterate over the partitions, which are odd numbers from 1 to
     geno_mask - 2 */
  for (unsigned int i = 1; i < geno_mask; i += 2)
    {
      /* iterate over the positions of the genome, testing to see where
	 there are changes */
      _Bool set_p = true;
      /* initialize the value to 1.0 */
      *result = 1.0;
      for (unsigned int j = 1; j < nloci; j++)
	{
	  /* save r's address, so we can iterate by modifying a pointer */
	  double * rptr = r;
	  /* if there is a change, we need to multiply */
	  if (bits_isset (i, j) != set_p)
	    *result *= *rptr;

	  /* otherwise we do nothing */
	  set_p = bits_isset (i, j);
	  /* advance our position in the r array */
	  rptr++;
	}
      /* now advance the array of results */
      result++;
    }
  return resultptr;
}

double
rec_total (size_t nloci, unsigned int j, unsigned int k,
	   unsigned int target, double * r)
{
  /* find the total probability of recombination over extended
     recombination map R, given parents J and K, and offspring
     TARGET */
  double total = 0.0;

  /* this is a mask for the whole genome, or the largest number
     represented by the genome */
  const unsigned int geno_mask = (1 << nloci) - 1;

  /* a pointer to reference elements of r */
  double * rptr = r;
  /* iterate over the partitions of the genome, represented by odd
     numbers from 1 to geno_mask - 2 */
  unsigned int comp = geno_mask - 1;
  for (unsigned int part = 1;
       part < geno_mask;
       /* increment partition PART by 2 (find the next odd number);
	  decrement complement COMP by 2 (find the next-lowest even
	  number)  */
       part += 2, comp -= 2, rptr++)
    {
      /* differences between genomes defined by partitions: in
	 Burger's notation I= part and J = comp, i is the target, j
	 and k are parents */
      unsigned int delta_iIjI = ((j & part) == (target & part))? 1 : 0;
      unsigned int delta_iIkI = ((k & part) == (target & part))? 1 : 0;
      unsigned int delta_iJkJ = ((k & comp) == (target & comp))? 1 : 0;
      unsigned int delta_iJjJ = ((j & comp) == (target & comp))? 1 : 0;

      /* does this recombination produce the target? */
      if ((delta_iIjI && delta_iJkJ)
	  || (delta_iIkI && delta_iJjJ))
	{
	  /* otherwise we need to find the probability of specific
	     recombinations that produce the target */
	  /* the recombination fraction for this partition */
	  double rI = *rptr;
	  total += (delta_iIjI * delta_iJkJ) * rI 
	    + (delta_iIkI * delta_iJjJ) * rI;
	}
      else continue;
    }
  /* we always need half the amount calculated, since there is another
     recombinant (or non-recombinant) offspring */
  return 0.5 * total;
}

rtable_t **
rec_gen_table (size_t nloci, size_t geno, double * r)
{
  /* first create rtable: an array of sparse matrices of length GENO */
  sparse_elt_t ** rtable = malloc (geno * sizeof (sparse_elt_t *));
  if (rtable == NULL)
    error (0, ENOMEM, "Null pointer\n");
  const unsigned int geno_mask = (1 << nloci) - 1;
  /* find the "extended recombination array" */
  double * xr = rec_extend_r (nloci, r);
  /* find the total probability of no recombination */
  double rnot = 1.0;
  for (int i = 0; i < (size_t) floor(geno_mask / 2); i++)
    rnot *= (1 - xr[i]);
  /* iterate over offspring entries, using endptr to keep track of
     position in the kth entry of rec_table, which is an array of
     GENO  */
  for (unsigned int target = 0; target< geno; target++)
    {   
      rtable[target] = sparse_new_elt (NULL, 0.0, NULL);
      sparse_elt_t * endptr = NULL;
      for (unsigned int k = 0; k < geno; k++)
	for (unsigned int j = 0; j < geno; j++)
	  {
	    _Bool new_elt_p = false;
	    double total;
	    unsigned int set;
	    unsigned int off;
	    unsigned int common;		
	    /* does the transpose already exist? */
	    if (isgreater(total = sparse_get_val (rtable[target], j, k),
			       0.0))
	      /* this avoids the function call to rec_total */
	      /* however, it doesn't work when the new entry *should*
		 be zero!  The later procedures should catch this, but
		 I'm fixing a bug where they don't...  */
	      new_elt_p = true;
	    else if ((target == j) && (j == k))
	      {
		/* do we need to make an entry? */
		new_elt_p = true;
		/* there is only one possibility: */
		total = 1.0;
	      }
	    else if (j == k)
	      /* no target offspring possible */
	      continue;
	    else if (target == k)
	      {
		/* Do one parent and the offspring have set bits in
		   common? */
		set = target & j;
		/* do they have off bits in common? */
		off = bits_extract(0, nloci, ~target & ~j);
		common = set | off;
		if (common)
		  /* if recombinant target offspring are possible,
		     then total = 0.5 */
		  total = 0.5;
		else if (common == 0)
		  total = 0.5 * rnot;
		else
		  total = rec_total (nloci, k, j, target, xr);
		if (total != 0.0)
		  new_elt_p = true;
		else
		  continue;
	      }
	    else if (target == j)
	      {
		/* Do one parent and the offspring have set bits in
		   common? */
		set = target & k;
		/* do they have off bits in common? */
		off = bits_extract(0, nloci, ~target & ~k);
		common = set | off;
		if (common)
		  /* if recombinant target offspring are possible,
		     then total = 0.5 */
		  total = 0.5;
		else if (common == 0)
		  total = 0.5 * rnot;
		else
		  /* if not, call rec_total */
		  total = rec_total (nloci, j, k, target, xr);
		if (total != 0.0)
		  new_elt_p = true;
		else
		  continue;
	      }
	    /* offspring must be recombinant or impossible */
	    else 
	      {
		/* Do one parent and the offspring have set bits in
		   common? */
		set = target & j;
		/* do they have off bits in common? */
		off = bits_extract(0, nloci, ~target & ~j);
		/* do they have any bits in common? */
		common = set | off;
		/* find where they are NOT alike: */
		unsigned int needed = bits_extract (0, nloci, ~common);
		/* if we need *all* the alleles from the other parent,
		   then we're screwed since we already know the other
		   parent is not identical to the target offspring */
		if (needed == geno_mask)
		  continue;
		/* now iterate over the needed bits; break at the
		   first allele we *can't* get from the other
		   parent */
		for (unsigned int pos = bits_ffs (needed);
		     /* remember that pos is 1-indexed! */
		     pos <= nloci && needed != 0;
		     pos = bits_ffs (needed))
		  {
		    /* do the bits from the other parent and the
		       offspring match? */
		    if (bits_isset (k, pos - 1) == bits_isset (target, pos - 1))
		      {
			/* for now the other parent is good */
			new_elt_p = true;
			/* clear that bit so we can reevaluate pos */
			needed &= ~(1 << (pos - 1));
		      }
		    else
		      {
			new_elt_p = false;
			break;
		      }
		  }
		if (new_elt_p)
		  total = rec_total (nloci, j, k, target, xr);
		else continue;
	      }	/* else */
	    
	    /* if we need a new element: */
	    if (new_elt_p && (endptr != NULL))
	      {
		/* create a link to the next element: */
		endptr->next = sparse_new_elt (NULL, total, NULL);
		/* increment endptr */
		endptr = endptr->next;	
	      }
	    /* if this is our first link in the chain: */
	    else if (new_elt_p)
	      {
		endptr = rtable[target];
		endptr->val = total;
	      }

	    /* set the indices */
	    endptr->indices[0] = k;
	    endptr->indices[1] = j;
	  }	/* for k < geno */
    }	    /* for target < geno */
  return rtable;
}

void
rec_mating (double * freqs, haploid_data_t * data)
{
  size_t geno = data->geno;
  sparse_elt_t ** rtable = data->rec_table;
  double ** mtable = data->mtable;
  /* find the frequencies of offspring from recombination table RTABLE
     and mating table MTABLE */

  /* FREQS[k] is the total of the Hadamard product of MTABLE and
     RTABLE[k] */
  for (int k = 0; k < geno; k++)
    freqs[k] = sparse_mat_tot (geno, mtable, rtable[k]);
}
