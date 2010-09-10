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
#include <stdio.h>
#include "sparse.h"
#include <float.h>
#include <assert.h>
#include <stdint.h>

#define PSWITCH(x) pholder = parent;		\
  parent = other;				\
  other = pholder;

double
rec_iterate (uint j, uint k, uint target, double * r, size_t nloci)
{
  /* degenerate condition (to stop recursion); recursion *should* stop
     when we reach the end of the while loop */
  if (nloci == 1) return 1.0;

  /* find the first locus where the parents are different from each
     other or the target */
  /* this is roach bait: check for bugs here first */
  uint parent, other;
  if (bits_isset (j, 0) == bits_isset (target, 0))
    {
      parent = j;
      other = k;
    }
  else
    {
      parent = k;
      other = j;
    }
  double total = 1.0F;
  do
    {
      /* go by pairs of loci */
      uint pmask = 3;
      uint tbits = target & pmask;
      uint pbits = parent & pmask;
      uint obits = other & pmask;
      /* how many bits do we need to change from the parents to make
	 the offspring? */
      int diff = bits_hamming (tbits, pbits) + bits_hamming (tbits, obits);
      uint pholder;		/* placeholder for switching parents */
      /* non-recombinant case */
      if ((diff == 2) && ((pbits == tbits)))
	total *= 1.0 - *r;
      /* recombinant case: */
      else if (diff == 2)
	{
	  total *= *r;
	  /* switch parents */
	  PSWITCH(void);
	}
      else if ((diff == 1) && (bits_hamming (tbits, pbits) == 1))
	{
	  /* if the difference arises from a difference from the
	     parent, then we should switch to the other parental
	     chromosome */
	  total *= *r;
	  /* switch parents */
	  PSWITCH(void);
	}
      else if ((diff == 1)
	       && (bits_isset (tbits, 1) != bits_isset (obits, 1)))
	/* then we're in a non-recombinant junction */
	total *= 1.0 - *r;
      else if (diff == 1)	/* branching point */
	{
	  /* recursion step */
	  total *= (*r * rec_iterate (other >> 1, parent >> 1, target >> 1,
				      r+1, nloci - 1)
		    + (1 - *r) * rec_iterate (parent >> 1, other >> 1, target >> 1,
					      r+1, nloci - 1));
	  break;
	}
      parent >>= 1; other >>= 1; target >>= 1; nloci--; 
    } while (nloci > 1);

   /* add up the expressions: */
  return total;
}

double
rec_total (uint j, uint k, uint target, double * r, size_t nloci)
{
  /* find the total probability of recombination over extended
     recombination map R, given parents J and K, and offspring
     TARGET */
  /* eliminate confounding cases right away */
  if ((j == k) && (j == target))
    /* only one possibility */
    return 1.0;
  else if ((j == k) || ((j ^ target) & (k ^ target)))
    /* no recombination is possible */
    return 0.0;
  else if ((bits_hamming (j, k) == 1) && ((j == target) || (k == target)))
    return 0.5;
  
  double result = rec_iterate (j, k, target, r, nloci); 
  if ((j & 1) == (k & 1))
    result += rec_iterate (k, j, target, r, nloci); 
 
  return result / 2.0;
}

void
rtable_new (rtable_t * rtable, double val, uint i, uint j)
{
  /* if the incoming rtable is NULL, then we need to create a fresh
     rtable object */
  if (rtable == NULL)
    rtable = sparse_new_elt (NULL, val, NULL);
  /* otherwise we need to create a new one and link it in to the old
     one, then advance it by one */
  rtable->val = val;
  rtable->indices[0] = i;
  rtable->indices[1] = j;
  /* set up the next link in the chain */
  int * indices = calloc (2, sizeof (int));
  if (indices == NULL)
    /* if you don't have enough space for two integers, we might as
       well give up: */
    error (0, ENOMEM, "Null pointer\n");
  rtable_t * new = sparse_new_elt (indices, 0, NULL);
  rtable->next = new;
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
  for (uint target = 0; target< geno; target++)
    {   
      rtable[target] = sparse_new_elt (NULL, 0.0, NULL);
      sparse_elt_t * endptr = rtable[target];
      for (uint k = 0; k < geno; k++)
	{
	  for (uint j = 0; j < geno; j++)
	    {
	      double total;
	      /* does the transpose already exist? */
	      if ((j == k) && (k == target))
		{
		  rtable_new (endptr, 1.0, k, j);
		  endptr = endptr->next;
		}
	      else if (isgreater(total = sparse_get_val (rtable[target], j, k), 0.0))
		{
		  rtable_new (endptr, total, k, j);
		  endptr = endptr->next;
		}
	      else if (isgreater(total = rec_total (k, j, target, r, nloci), 0.0))
		{
		  rtable_new (endptr, total, k, j);
		  endptr = endptr->next;
		}
	      else continue;
	    }
	}      /* for k < geno */
    } /* for target < geno */
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
