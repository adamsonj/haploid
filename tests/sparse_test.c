#include "../src/haploidtest.h"


#define LEN 0xf

int
main (void)
{
  /* matrix to multiply times identity matrix */
  double ** matelts = calloc (LEN, sizeof (double *));
  if (matelts == NULL) error (0, ENOMEM, "Null pointer\n");
  
  /* identity matrix: */
  int matindx[LEN][2];
  /* array of matrix elements */
  sparse_elt_t * mat[LEN];
  int i, j;

  for (i = 0; i < LEN; i++)
    {
      for (j = 0; j < 2; j++)
	matindx[i][j] = i;

      mat[i] = sparse_new_elt (matindx[i], 1.0, NULL);
      /* set the ->next of each element, except the last one, which
	 should stay NULL */
      if (i > 0)
	mat[i-1]->next = mat[i];

      /* now fill in the matrix for multiplication */
      matelts[i] = calloc (LEN, sizeof (double));
      if (matelts[i] == NULL)
	error (0, ENOMEM, "Null pointer\n");

      for (j = 0; j < LEN; j++)
	matelts[i][j] = i;
    }

  sparse_elt_t * matresult = sparse_mat_mat_kron (LEN, matelts, *mat);
  for (sparse_elt_t * matptr = matresult;
       matptr->next != NULL;
       matptr = matptr->next)
    printf ("M[%2i,%2i] = %9.8f\n",
	    /* row index */
	    matptr->indices[0],
	    /* column index */
	    matptr->indices[1],
	    /* value */
	    matptr->val);

 
  double res;
  res = sparse_mat_tot (matresult);

  printf ("Total of entries in M = %9.8f \n", res);

  return 0;
}
	 
