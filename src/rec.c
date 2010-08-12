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
  const uint geno_mask = (1 << nloci) - 1;
  double * result = malloc ((size_t) (geno_mask >> 1) * sizeof (double));
  if (result == NULL)
    error (0, ENOMEM, "Null pointer\n");
  /* we will access the values of result by pointer arithmetic (to
     avoid calculating indices) */
  double * resultptr = result;

  /* iterate over the partitions, which are odd numbers from 1 to
     geno_mask - 2 */
  for (uint i = 1; i < geno_mask; i += 2)
    {
      /* iterate over the positions of the genome, testing to see where
	 there are changes */
      _Bool set_p = true;
      /* initialize the value to 1.0 */
      *result = 1.0;
      for (uint j = 1; j < nloci; j++)
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
rec_total (size_t nloci, uint j, uint k,
	   uint target, double * r)
{
  /* find the total probability of recombination over extended
     recombination map R, given parents J and K, and offspring
     TARGET */
  double total = 0.0;
  const uint geno_mask = (1 << nloci) - 1;
  double * rptr = r;
  /* iterate over the partitions of the genome, represented by odd
     numbers from 1 to geno_mask - 2 */
  uint comp = geno_mask - 1;
  /* increment partition PART by 2 (find the next odd number) */
  for (uint part = 1; part < geno_mask; part += 2,  comp -= 2, rptr++)
    {
      /* differences between genomes defined by partitions: in
	 Burger's notation I= part and J = comp, i is the target, j
	 and k are parents */
      uint delta_iIjI = ((j & part) == (target & part))? 1 : 0;
      uint delta_iIkI = ((k & part) == (target & part))? 1 : 0;
      uint delta_iJkJ = ((k & comp) == (target & comp))? 1 : 0;
      uint delta_iJjJ = ((j & comp) == (target & comp))? 1 : 0;

      /* does this recombination produce the target? */
      if ((delta_iIjI && delta_iJkJ) || (delta_iIkI && delta_iJjJ))
	total += *rptr;
      else continue;
    }
  return 0.5 * total;
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

int
rec_case (uint target, uint mom, uint dad,
	  rtable_t * rtable, double * r, uint nloci)
{
  const uint	genomask = (1 << nloci) - 1;
  uint	setmom	 = target & mom; uint	setdad	 = target & dad;
  uint	offmom	 = (~target & genomask) & (~mom & genomask);
  uint	offdad	 = (~target & genomask) & (~dad & genomask);
  /* is recombination possible? */
  /* common alleles between mom and target */
  uint	momcom	 = setmom | offmom; uint dadcom	 = setdad | offdad; 
  /* alleles needed by mom to meet target */
  uint	needmom	 = ~momcom & genomask; uint needdad = ~dadcom & genomask; 
  /* if non-zero, recombination with dad genome can yield target */
  int momp = needmom & dadcom; int dadp = needdad & momcom; 
  /* find the "extended recombination array" */
  double * xr = rec_extend_r (nloci, r);
  if ((mom == dad) && (dad == target))
    rtable_new (rtable, 1.0, mom, dad);
  else if (mom == dad)
    /* no new element */
    return 1;
  else if ((!needmom && !dadcom) || (!needdad && !momcom))
    /* mom (dad) equals target and dad (mom) has no genes to offer;
       therefore we need the probability of no recombination */
    {
      double rnot = 1.0;
      for (int i = 0; i < (genomask >> 1) ; i++)
	rnot *= (1 - xr[i]);	/* no recombination */
      rtable_new (rtable, 0.5 * rnot, mom, dad);
    }
  else if ((needmom && !dadcom) || (needdad && !momcom))
    /* no recombination possible */
    return 1;
  else if (((needmom | dadp) == dadp) || ((needdad | momp) == momp))
    /* mom (dad) = target, but can recombine with dad (mom) to produce
       recombinant target */
    rtable_new (rtable, 0.5, mom, dad);
  else if (momp)
    /* get the total recombination value and create a new element */
    {
      double total = rec_total (nloci, mom, dad, target, xr);
      rtable_new (rtable, total, mom, dad);
    }      
  else
    return 1;  
  return 0;
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
	    if (isgreater(total = sparse_get_val (rtable[target], j, k), 0.0))
	      {
		rtable_new (endptr, total, k, j);
		endptr = endptr->next;
	      }
	    else
	      {
		int rc_ck = rec_case (target, k, j, endptr, r, nloci);
		if (!rc_ck)
		  endptr = endptr->next;
		else continue;
	      }		
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
