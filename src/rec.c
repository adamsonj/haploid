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
#include "haploidpriv.h"

double
rec_total (unsigned int nloci, unsigned int diff,
	   double * r, _Bool recomb_p)
{
  /* find the total probability of recombination (or not, depending on
     RECOMB_P) over recombination map R, given difference map DIFF */
  double total = 1.0;
  double * rptr = r;		/* iterate this pointer */
  /* skip ahead to the first set bit: */
  unsigned int pos = bits_ffs (diff);
  unsigned int subpos;
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

  return total;
}

rtable_t **
rec_gen_table (int nloci, int geno, double * r)
{
  /* first create rtable: an array of sparse matrices of length GENO */
  /* using malloc () here enables us to put rtable on the heap, so
     that we can return its address, thereby returning a pointer to a
     sparse matrix; this allows for a more functional programming
     style; i.e. creates the illusion of generating an array of sparse
     matrices for the user */
  sparse_elt_t ** rtable = malloc (geno * sizeof (sparse_elt_t *));
  if (rtable == NULL)
    error (0, ENOMEM, "Null pointer\n");
  unsigned int i;

  double sum = 0;
  for (i = 0; i < nloci - 1; i++)
    sum += 1 - r[i];
  /* iterate over offspring entries, using endptr to keep track of
     position in the kth entry of rec_table, which is an array of
     GENO  */
  sparse_elt_t * endptr;
  for (unsigned int k = 0; k < geno; k++)
    {   
      rtable[k] = sparse_new_elt (NULL, 0.0, NULL);
      endptr = rtable[k];
      for (i = 0; i < geno; i++)
	for (unsigned int j = 0; j < geno; j++)
	  {
	    if ((i == j) && (j == k))
	      /* there is only one possibility: */
	      endptr->val = 1.0;
	    else if (i == j)
	      /* no target offspring possible */
	      continue;
	    else if (j == k)
	      {
		/* this case is (1-r)/2 */
		endptr->val = 0.5 * rec_total (nloci, i ^ k, r, false);
		/* are recombinant offspring possible? i.e. is the
		   other parent different at all sites? */
		if (i != bits_extract (0, nloci, ~k)) 
		  endptr->val += 0.5 * rec_total (nloci, i ^ k, r, true);
	      }
	    else if (i == k)
	      {
		/* this case is (1-r)/2 */
		endptr->val = 0.5 * rec_total (nloci, j ^ k, r, false);
		if (j != bits_extract (0, nloci, ~k))
		  /* recombinant offspring are also possible */
		  endptr->val += 0.5 * rec_total (nloci, j ^ k, r, true);
	      }
	    else if ((i != bits_extract (0, nloci, ~k))
		     && (j != bits_extract (0, nloci, ~k)))
	      /* offspring must be recombinant */
	      endptr->val = 0.5 * rec_total (nloci, i ^ k, r, true);
	    else
	      /* if they have no alleles in common there is no need to
		 set the indices or advance endptr; proceed to the top
		 of the innermost loop */
	      continue;

	    /* we must have some alleles in common; if not we do not
	       make an entry in the table */
	    endptr->indices[0] = i;
	    endptr->indices[1] = j;
	    
	    /* create a link to the next element: */
	    endptr->next = sparse_new_elt (NULL, 0.0, NULL);
	    /* increment endptr */
	    endptr = endptr->next;
	  }
    }
  return rtable;
}

double *
rec_mating (haploid_data_t * data)
{
  int geno = data->geno;
  sparse_elt_t ** rtable = data->rec_table;
  double ** mtable = data->mtable;
 
  /* find the frequencies of offspring from recombination table RTABLE
     and mating table MTABLE */

  /* a vector for results: */
  double * result = calloc (geno, sizeof (double));
  /* RESULT[k] is the total of the Kronecker product of MTABLE and
     RTABLE[k] */
  int k;
  sparse_elt_t * recmat[geno];
  /* first an OpenMP directive: */
#pragma omp parallel for shared(rtable, mtable, recmat, result)
  for (k = 0; k < geno; k++)
    {
      recmat[k] = sparse_mat_mat_kron (geno, mtable, rtable[k]);
      result[k] = sparse_mat_tot (recmat[k]);
    }

  return result;
}
