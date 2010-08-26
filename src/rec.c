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

uint
rec_more_zygs (uint stem, uint nloci)
{
  /* iterate over the integers (genotypes), &-ing them with stem,
     returning the total number of zygotes with the stem */
  int n = 0;
  /* set up the stem */
  uint genomask = (1 << nloci) - 1;
  uint mask = ~stem & genomask;
  for (int i = 0; i < (1 << nloci); i++, mask &= i) 
    if ((mask & i) != mask) n++;
  return n;
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

  /* number of recombination sites */
  uint njunx = nloci - 1;
  double total = 1.0;

  uint p = 0; 
  for (; p < njunx - 1; p++)
    {
      /* go by pairs of loci */
      uint jmask = 0x3 << p;
      uint mt = target & jmask;
      uint jt = j & jmask;
      uint kt = k & jmask;
      int diff = bits_hamming (jt, mt) + bits_hamming (kt, mt);

      /* non-recombinant case */
      if ((diff == 2) && ((jt == mt) || (kt == mt)))
	total *= (1 - r[p]);
      /* recombinant case: */
      else if (diff == 2)
	total *= r[p];
    }
  return total / 2.0;
}

void
rtable_new (rtable_t * rtable, double val, uint i, uint j)
{
  /* if the incoming rtable is NULL, then we need to create a fresh
     rtable object */
  if (rtable == NULL)
    {
      rtable = malloc (sizeof (rtable_t));
      if (rtable == NULL)
	error (0, ENOMEM, "Null pointer\n");	
    }
  /* otherwise we need to create a new one and link it in to the old
     one, then advance it by one */
  rtable->val = val;
  rtable->indices[0] = i;
  rtable->indices[1] = j;
  /* set up the next link in the chain */
  rtable_t * new = malloc (sizeof (rtable_t));
  new->indices = calloc (2, sizeof (uint));
  if ((new == NULL) || (new->indices == NULL))
    error (0, ENOMEM, "Null pointer\n");
  else
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
	      if (j == k && k == target)
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
