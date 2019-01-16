#include <stdio.h>
#include "obsidian.h"

int
print_array(double * array, int length) {
  for (int index = 0; index < length; index++)
    printf("array[%d] = %f\n", index, array[index]);
  return 0;
}

evolve_result
evolve(int reactions_length,
       int substrates_length,
       double * stoichiometry,
       double * rates,
       double * state,
       double duration) {
  evolve_result result = {NULL, NULL, NULL};
  return result;
}
