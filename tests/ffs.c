#include "../src/haploidtest.h"

int
main (void)
{
  for (int i = 0; i < 0x11; i++)
    printf ("bits_ffs(%2x) = %4x\n", i, bits_ffs (i));
  return 0;
}
  
