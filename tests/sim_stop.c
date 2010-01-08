#define GENO 4

#include "../haplib/haplib.h"

#define LEN 2

int
main (void)
{
  /*

    produce a small, silly simulation that will stop after a few
    rounds

  */

  
  /* define an array */
  double simarr[LEN] = {1.0, 0.5};
  double oldarr[LEN] = {0.0, 0.0};  
  int n = 2;

  /* counter for loop */
  int i;

  while (sim_stop_ck (simarr, oldarr, LEN, 1e-16) && n < 1e6)
    {

      /*

	while you reduce the size of both elements to zero, check them
	successively using sim_stop_ck; use large changes at first, then
	reduce them to zero

      */
    
    
      printf ("x(%2d) = [%16.15f, %16.15f]\n", n-1, simarr[0], simarr[1]);
      for (i = 0; i < LEN; i++)
	{
	  oldarr[i] = simarr[i];

	  simarr[i] = simarr[i] / (pow (n, 2));
	}

      n++;
    
    };
  if (n < 1e6)
    return 0;
  else
    return 1;
}
