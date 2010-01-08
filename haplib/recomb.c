/* declarations */
#include "haplibpriv.h"


/* note: this should be replaced with a bit-hacking macro */
static inline int
get_bit (int integer, int i)
{
  /*  get_bit returns the ith bit of integer, i = 0,1,... */
  int pw21 = (int) pow (2, i+1);
  int pw2 = (int) pow (2,i);
  return (integer % pw21) / pw2;
}

static double
rec_prob (int mask, int nloci, double * r)
{
  /*  rec_prob calculates the probability of a given recombination pattern. */
  /*  The pairwise recombination rates r[0],r[1],...,r[nLoci-1] are
      declared and initialized external to this routine.  */
  /*  M.K., II-96.  */
  /*  DEBUGGED */
  double prob = 1.0;
  int i;

  for (i = 0; i < nloci - 1; i++)
    {
      if (get_bit (mask, i) == get_bit (mask, i + 1))
	prob *= (1.0 - r[i]);
      else
	prob *= r[i];
    }
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

  unsigned int zyg1 = 0, zyg2 = 0;
  unsigned int bit[nloci];
  int k;
  int powk;

  /* preliminaries         */
  for (k = 0; k < nloci; k++)
    {

      powk = pow (2, k);

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
  *zyg = zyg1;
  zyg++;
  *zyg = zyg2;
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
  const int mask_lim = (int) pow (2, nloci - 1);
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

	  /* if j != i, then we have two cases */

	  /* if j < i then we can mirror the above-diagonal entries */
	  if (j < i)
	    {
	    for (k = 0; k < geno; k++)
	      rec_table[i][j][k] = rec_table[j][i][k];
	    }
	  /* if j > i then we need to calculate the entry; oh darn! */
	  else
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
		  /* if there are not enough entries in r, zygote_prob
		     could return nan */
		  if (isunordered (zygote_prob, 0.0))
		    {
		      perror ("zygote_prob is not a number\n");
		      return FE_INVALID;
		    }
		  else
		    {
		    

		      /* Generate zygote genotypes
		 
			 We should get two out of using any particular mask:
			 the one from specified directly by mask, and the one
			 specified by its its bitwise complement; store the
			 result in zyg  */

		      zygote_genotypes (i, j, mask, zyg, geno);

		      /* Find indices corresponding to the zygotes produced by
			 mask and parents, and add to their probabilities, as
			 we just found a new way to generate them  */

		      for (k = 0; k < 2; k++)
			rec_table[i][j][zyg[k]] += zygote_prob;
		    }
		}
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
