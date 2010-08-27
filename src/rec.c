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
  uint H = bits_hamming (j, k);
  /* eliminate confounding cases right away */
  if ((j == k) && (j == target))
    /* only one possibility */
    return 1.0;
  else if ((j == k) || ((j ^ target) & (k ^ target)))
    /* no recombination is possible */
    return 0.0;
  else if ((H == 1) && ((j == target) || (k == target)))
    return 0.5;

  /* number of recombination sites */
  uint njunx = nloci - 1;
  size_t neqns = 1 << njunx;
  double total[neqns];
  total[0] = 1.0;
  /* the spot to insert a new entry */
  size_t slot = 1;
  for (int i = slot; i < neqns; i++) total[i] = 0.0;
  size_t maybes = 0;
  uint p = 0; 
  do
    {
      /* go by pairs of loci */
      uint jmask = 0x3 << p;
      uint mt = target & jmask;
      uint jt = j & jmask;
      uint kt = k & jmask;
      /* how many bits do we need to change from the parents to make
	 the offspring? */
      int diff = bits_hamming (jt, mt) + bits_hamming (kt, mt);

      /* non-recombinant case */
      if ((diff == 2) && ((jt == mt) || (kt == mt)))
	for (int i = 0; i < slot; i++) total[i] *= (1 - r[p]);
      /* recombinant case: */
      else if (diff == 2)
	for (int i = 0; i < slot; i++) total[i] *= r[p];
      else if (diff == 1)
	{
	  /* expand the array; new entries (slot to slot*2) are old
	     entries multiplied by r, then multiply old entries (the
	     first slot of them) by 1-r */
	  for (int i = 0; i < slot; i++)
	    {
	      total[slot + i] = total[i] * r[p];
	      total[i] *= 1 - r[p];
	    }
	  /* we may need to check for this value going over neqns */
	  slot *= 2;
	  maybes++;
	}
      p++;
    } while (p < njunx); 

  if ((j != target) && (k != target) && (maybes == njunx))
    /* offspring must be recombinant, so we need to eliminate the
       possibility of no recombination when every junction is a
       maybe */
    total[0] = total[neqns - 1] = 0;
  else if (maybes == njunx)
    /* OTOH if one of the parents does equal the target and we have
       all maybes then we want only the sum of the totally recombinant
       or totally nonrecombinant (which do not add to 1
       incidentally!) */
    for (int i = 1; i < neqns - 1; i++) total[i] = 0.0;

  /* add up the expressions: */
  double result = 0.0;
  for (int i = 0; i < neqns; i++) result += total[i];    
  return result / 2.0;
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
