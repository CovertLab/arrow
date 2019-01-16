#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "obsidian.h"

int
print_array(double * array, int length) {
  for (int index = 0; index < length; index++)
    printf("array[%d] = %f\n", index, array[index]);
  return 0;
}

static int INITIAL_LENGTH = 1000;

int
choose(int n, int k) {
  double combinations = 1.0;
  int i;
  for (i = 0; i < k; i++) {
    combinations *= (n - i) / (i + 1);
  }
  return combinations;
}

evolve_result
evolve(int reactions_length,
       int substrates_length,
       double * stoichiometry,
       double * rates,

       int * reactants_lengths,
       int * reactants_indexes,
       int * reactants,
       double * reactions,
       
       int * dependencies_lengths,
       int * dependencies_indexes,
       int * dependencies,

       double duration,
       double * state) {

  double * time = malloc((sizeof (double *)) * INITIAL_LENGTH);
  int * events = malloc((sizeof (int *)) * INITIAL_LENGTH);

  double * propensities = malloc((sizeof (double *)) * reactions_length);
  int * update = malloc((sizeof (int *)) * reactions_length);
  int update_length = reactions_length;

  int reaction, reactant, species, index;
  double total, interval, sample, progress, point;
  
  int choice, step = 0, up = 0;
  double now = 0.0;

  for (reaction = 0; reaction < reactions_length; reaction++) {
    update[reaction] = reaction;
  }

  while(now < duration) {
    printf("step %i", step);

    for (up = 0; up < update_length; up++) {
      int reaction = update[up];
      propensities[reaction] = rates[reaction];

      for (reactant = 0; reactant < reactants_lengths[reaction]; reactant++) {
        index = reactants_indexes[reaction] + reactant;
        propensities[reaction] *= choose(state[reactants[index]], reactions[index]);
      }
    }

    printf("propensities");
    print_array(propensities, reactions_length);

    total = 0.0;
    for (reaction = 0; reaction < reactions_length; reaction++) {
      total += propensities[reaction];
    }

    printf("total: %f\n", total);

    if (total == 0) {
      interval = 0.0;
      choice = -1;
    } else {
      sample = (double) rand() / RAND_MAX;
      interval = -log(1 - sample) / total;
      point = ((double) rand() / RAND_MAX) * total;

      printf("sample: %f - interval: %f - point: %f\n", sample, interval, point);

      choice = 0;
      progress = 0.0;
      while(progress + propensities[choice] < point) {
        progress += propensities[choice];
        choice += 1;
      }

      printf("progress: %f - choice: %d - now: %f - duration: %f\n", progress, choice, now, duration);

      if ((choice == -1) || ((now + interval) > duration)) {
        printf("completing: %d %d", choice == -1, (now + interval) > duration);
        break;
      }

      now += interval;

      printf("new time: %f\n", now);

      time[step] = now;
      events[step] = choice;

      for (species = 0; species < substrates_length; species++) {
        state[species] += stoichiometry[choice * substrates_length + species];
      }

      update_length = dependencies_lengths[choice];
      for (up = 0; up < update_length; up++) {
        update[up] = dependencies[dependencies_indexes[choice] + up];
      }

      printf("update\n");
      print_array(update, update_length);
    }

    step += 1;
  }

  print_array(state, substrates_length);
  
  evolve_result result = {0, time, events, state};

  free(propensities);
  free(update);

  return result;
}
