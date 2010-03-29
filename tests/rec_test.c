#include "../src/haploidtest.h"
#include <time.h>

#define NLOCI 4
#define GENO 16
#define LEN 3
double r[LEN] = { 0.25, 0.125, 0.0625};

rtable_t ** rec_table;

/* declarations: */
void
rec_test_total (void);

void
rec_test_prtable (rtable_t ** rtable);

int
main (void)
{
  rec_test_total ();
  rec_table = rec_gen_table (NLOCI, GENO, r);
  rec_test_prtable (rec_table);
  double freq[GENO];
  /* time the next calculation */
  for (int j = 0; j < GENO; j++)
    freq[j] = 1.0 / GENO;
  haploid_data_t rec_test_data =
    { GENO, NLOCI, rec_table, rmtable (GENO, freq)};
  time_t time1, time2;
  time1 = time (NULL);
  rec_mating (&rec_test_data);
  time2 = time (NULL);
  /* print freq */
  for (int j = 0; j < GENO; j++)
    printf ("x[%1x] = %9.8f\n", j, freq[j]);

  printf ("Call to rec_mating () took %9.8f sec\n", difftime(time1, time2));
  return 0;
}

void
rec_test_total (void)
{
  printf ("Recombination probabilities: ");
  for (int i = 0; i < LEN; i++)
    printf ("%9.8f ", r[i]);
  printf ("\n");
  printf ("Probability of recombinant offspring 0x1 from parent 0xF: %9.8f\n",
	  0.5 * rec_total (NLOCI, (0x0 ^ 0x1), r, true));
  printf ("Probability of non-recombinant offspring 0x1 \n"
	  "from parents 0xF X 0x1: %9.8f\n",
	  0.5 * rec_total (NLOCI, (0xF ^ 0x1), r, false));
  printf ("Rock 'n roll dynamite, baby!.\n");
}

/* print a recombination table */

void
rec_test_prtable (rtable_t ** rtable)
{
  /* traverse the table, printing zeros where there are no entries */
  sparse_elt_t * tptr;
  int geno;			/* offspring genotype */
  double val;
  int i, j;			/* row, column indices */
  for (geno = 0; geno < GENO; geno++)
    /* iterate over the rows of the table, printing the whole thing as
       a matrix (with zero entries) */
    {
      tptr = rtable[geno];

      /* print a header announcing the genotype: */
      printf ("Offspring %x:\n", geno);
      /* print the parent genotypes across the top */
      for (j = -1; j < GENO; j++)
	if (j < 0)
	  printf ("%2s", " ");
	else
	  printf ("%10x ", j);
      printf ("\n");
      for (i = 0; i < GENO; i++)
	{
	  /* print the paternal genotype for the row */
	  printf ("%x:", i);
	  /* now print probabilities */
	  for (j = 0; j < GENO; j++)
	    {
	      val = sparse_get_val (tptr, i, j);
	      printf ("%9.8f ", val);
	    }
	  printf ("\n");
	}
      printf ("\n");
    }	  
}
