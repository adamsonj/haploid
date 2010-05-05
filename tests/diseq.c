/*

  diseq.c: testing linkage disequilibrium function 

  Copyright 2010 Joel J. Adamson

  $Id$

  Joel J. Adamson -- http://www.unc.edu/~adamsonj
  University of North Carolina at Chapel Hill
  CB #3280, Coker Hall
  Chapel Hill, NC 27599-3280 <adamsonj@email.unc.edu>

  This file is part of haploid

  haploid is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  haploid is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with haploid.  If not, see <http://www.gnu.org/licenses/>.

*/

/* Commentary:

 This file tests the function ld_from_geno ().  The method by hand is
 either the typical Felsenstein linkagedisequilibrium x1*x2 - x3*x4 or
 a mimic of the routine employed by the function; a better test might
 be to throw random numbers into the function to see how often it
 produces a clearly wrong value (e.g. greater than 1).

 If you define the number of loci using -DNLOCI, make sure you define
 the proper number of genotypes using -DGENO.  Otherwise this will run
 with the canned values of 4 and 16.  I recommend testing this with
 TLTA values since they are easy to check literally "by hand."
  
*/

#include "../src/haploidtest.h"
#ifndef NLOCI
#define NLOCI 4
#endif

#ifndef GENO
#define GENO 16
#endif

int
main (void)
{
  const size_t nloci = NLOCI;
  const size_t genotypes = GENO;
  const double geno_zero =  M_PI_4;
  double allele_freqs[NLOCI];
  double genotype_freqs[GENO] = { [0] = geno_zero };
  for (int i = 0; i < genotypes; i++)
  {
    if (i == 0);
    else genotype_freqs[i] = (1.0 - geno_zero) / (genotypes - 1);
#ifdef DEBUG
    fprintf (stdout, "x[%2x] = %f\n", i, genotype_freqs[i]);
#endif
  }
  
  genotype_to_allele (allele_freqs, genotype_freqs, nloci, genotypes);
#ifdef DEBUG
  fprintf (stdout, "Allele frequencies:\n");
  
  for (int k = 0; k < genotypes; k++)
      fprintf (stdout, "x[%2x] = %f\n", k, allele_freqs[k]);
#endif
  /* You can't break up with me, I've got hand! */
  double hand;
  
#if (NLOCI > 2)
  double prod = 1.0;
  for (int j = 0; j < nloci; j++)
    prod *= 1 - allele_freqs[j];
  /* calculate the value explicitly */
  hand = geno_zero - prod;
#elif (NLOCI == 2)
  hand = genotype_freqs[0] * genotype_freqs[3]
    - genotype_freqs[1] * genotype_freqs[2];
#endif

  /* calculate the value with the library routine: */
  double ld = ld_from_geno (nloci, genotypes, genotype_freqs);
  fprintf (stdout, "Value by hand: %f\n", hand);
  fprintf (stdout, "Value by ld_from_geno: %f\n", ld);
  _Bool nomatch = (int)(hand - ld);
  if (nomatch) exit (EXIT_FAILURE);
  exit (EXIT_SUCCESS);
}
/* end of diseq.c */
