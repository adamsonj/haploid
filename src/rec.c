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

double
rec_total (size_t nloci, unsigned int diff,
	   double * r, _Bool recomb_p)
{
  /* find the total probability of recombination (or not, depending on
     RECOMB_P) over recombination map R, given difference map DIFF */
  /* WARNING: this function assumes that the other parent (the one NOT
     used to find DIFF provides the right alleles; the caller must
     verify this assumption) */
  double total = 1.0;
  double * rptr = r;		/* iterate this pointer */
  /* skip ahead to the first set bit: */
  unsigned int pos = bits_ffs (diff);
  unsigned int subpos;
  if (diff == 0)		/* special case */
    for (pos = 1; pos < nloci; rptr++, pos++)
      total *= recomb_p?(*rptr):(1 - *rptr);
  else
    {
      if (pos == 1)			/* special case */
	{
	  diff = ~diff;
	  pos = bits_ffs (diff);
	}
      do
	{
	  rptr += pos - 2;
	  total *= recomb_p?(*rptr):(1 - *rptr);
	  /* we know the bit at index pos-1 is set, so if the bit at index
	     pos is set, there is no transition and we can skip ahead; or
	     we can find the next unset bit by usign bits_ffs on the one's
	     complement! */
	  subpos = bits_ffs(bits_extract (pos, nloci, ~diff));
	  /* if subpos is 0, then there are no more transitions and our
	     work here is done */
	  if (subpos == 0)
	    break;
	  else
	    /* find the next transition */
	    pos += subpos; 
	}
      while (pos < nloci);
    }
  return total;
}

rtable_t **
rec_gen_table (size_t nloci, size_t geno, double * r)
{
  /* first create rtable: an array of sparse matrices of length GENO */
  sparse_elt_t ** rtable = malloc (geno * sizeof (sparse_elt_t *));
  if (rtable == NULL)
    error (0, ENOMEM, "Null pointer\n");

  /* iterate over offspring entries, using endptr to keep track of
     position in the kth entry of rec_table, which is an array of
     GENO  */
  for (unsigned int k = 0; k < geno; k++)
    {   
      rtable[k] = sparse_new_elt (NULL, 0.0, NULL);
      sparse_elt_t * endptr = NULL;
      for (unsigned int i = 0; i < geno; i++)
	for (unsigned int j = 0; j < geno; j++)
	  {
	    _Bool new_elt_p = false;
	    double total;
	    /* does the transpose already exist? */
	    if ((total = sparse_get_val (rtable[k], j, i)) != 0.0)
	      /* this avoids the function call to rec_total */
	      /* however, it doesn't work when the new entry *should*
		 be zero!  The later procedures should catch this, but
		 right I'm fixing a bug where they don't...  */
	      new_elt_p = true;
	    else if ((i == j) && (j == k))
	      {
		/* do we need to make an entry? */
		new_elt_p = true;
		/* there is only one possibility: */
		total = 1.0;
	      }
	    else if (i == j)
	      /* no target offspring possible */
	      continue;
	    else if (j == k)
	      {
		new_elt_p = true;
		/* this case is (1-r)/2 */
		total = 0.5 * rec_total (nloci, 0, r, false)
		  /* recombinant offspring are also possible: we know
		     we can get the alleles we need from the parent
		     identical to the offspring */
		  + 0.5 * rec_total (nloci, i ^ k, r, true);
	      }
	    else if (i == k)
	      {
		new_elt_p = true;
		/* this case is (1-r)/2 */
		total = 0.5 * rec_total (nloci, 0, r, false)	
		  /* recombinant offspring are also possible: we know
		     we can get the alleles we need from the parent
		     identical to the offspring */
		  + 0.5 * rec_total (nloci, j ^ k, r, true);
	      }
	    /* offspring must be recombinant or impossible */
	    else 
	      {
		/* Do the parent and the offspring have set bits
		   in common? */
		unsigned int set = i & k;
		/* do they have off bits in common? */
		unsigned int off = bits_extract(0, nloci, ~i & ~k);
		/* do they have any bits in common? */
		unsigned int common = set | off;
		/* find where they are NOT alike: */
		unsigned int needed = bits_extract (0, nloci, ~common);
		/* if we need *all* the alleles from the other parent,
		   then we're screwed since we already know the other
		   parent is not identical to the target offspring */
		if (needed == (1 << nloci) - 1)
		  continue;
		/* now iterate over the needed bits; break at the
		   first allele we *can't* get from the other
		   parent */
		for (unsigned int pos = bits_ffs (needed);
		     pos < nloci && needed != 0;
		       pos = bits_ffs (needed))
		  {
		    /* do the bits from the other parent and the
		       offspring match? */
		    if (bits_isset (j, pos - 1) == bits_isset (k, pos - 1))
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
		  total = 0.5 * rec_total (nloci, i ^ k, r, true);
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
		endptr = rtable[k];
		endptr->val = total;
	      }

	    /* set the indices */
	    endptr->indices[0] = i;
	    endptr->indices[1] = j;
	  }	/* for j < geno */
    }	    /* for k < geno */
  return rtable;
}

void
rec_mating (double * freqs, haploid_data_t * data)
{
  size_t geno = data->geno;
  sparse_elt_t ** rtable = data->rec_table;
  double ** mtable = data->mtable;
  double denom = 0.0;
  /* find the frequencies of offspring from recombination table RTABLE
     and mating table MTABLE */

  /* FREQS[k] is the total of the Kronecker product of MTABLE and
     RTABLE[k] */
  for (int k = 0; k < geno; k++)
    {
      freqs[k] = sparse_mat_tot (geno, mtable, rtable[k]);
      denom += freqs[k];
    }
  /* to ensure that frequencies sum to 1 */
  for (int k = 0; k < geno; k++)
    freqs[k] /= denom;
}
