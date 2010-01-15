#include "../src/haploid.h"

int
main (void)
{
  int i;
  int lim;
  lim = (int)pow (2,sizeof(unsigned int));
  for (i = 0; i < lim; i++)
    printf ("There are %d one-bits in %d.\n", POP(i), i);
  return 0;
}
