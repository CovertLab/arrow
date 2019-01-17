#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "obsidian.h"

int
print_array(double * array, int length) {
  for (int index = 0; index < length; index++) {
    printf("a[%d] = %f", index, array[index]);
    if (index == length - 1) {
      printf("\n");
    } else {
      printf(", ");
    }
  }

  return 0;
}

int
print_long_array(long * array, int length) {
  for (int index = 0; index < length; index++) {
    printf("a[%d] = %ld", index, array[index]);
    if (index == length - 1) {
      printf("\n");
    } else {
      printf(", ");
    }
  }

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

       long * reactants_lengths,
       long * reactants_indexes,
       long * reactants,
       double * reactions,
       
       long * dependencies_lengths,
       long * dependencies_indexes,
       long * dependencies,

       long * actors_lengths,
       long * actors_indexes,
       long * actors,

       double duration,
       double * state) {

  double * time = malloc((sizeof (double *)) * INITIAL_LENGTH);
  long * events = malloc((sizeof (long *)) * INITIAL_LENGTH);
  double * outcome = malloc((sizeof (double *)) * substrates_length);

  double * propensities = malloc((sizeof (double *)) * reactions_length);
  long * update = malloc((sizeof (long *)) * reactions_length);
  long update_length = reactions_length;
  long actors_length;

  long reaction, reactant, species, index, actor;
  double total, interval, sample, progress, point, adjustment, count;
  
  int choice, step = 0, up = 0;
  double now = 0.0;

  for (species = 0; species < substrates_length; species++) {
    outcome[species] = state[species];
  }

  for (reaction = 0; reaction < reactions_length; reaction++) {
    update[reaction] = reaction;
  }

  /* int rlength = 0; */
  /* for (reaction = 0; reaction < reactions_length; reaction++) { */
  /*   rlength += reactants_lengths[reaction]; */
  /* } */

  /* int dlength = 0; */
  /* for (reaction = 0; reaction < reactions_length; reaction++) { */
  /*   dlength += dependencies_lengths[reaction]; */
  /* } */

  /* printf("reactants_lengths: "); */
  /* print_long_array(reactants_lengths, reactions_length); */
  /* printf("reactants_indexes: "); */
  /* print_long_array(reactants_indexes, reactions_length); */
  /* printf("reactants: "); */
  /* print_long_array(reactants, rlength); */
  /* printf("reactions: "); */
  /* print_array(reactions, rlength); */

  /* printf("dependencies_lengths: "); */
  /* print_long_array(dependencies_lengths, reactions_length); */
  /* printf("dependencies_indexes: "); */
  /* print_long_array(dependencies_indexes, reactions_length); */
  /* printf("dependencies: "); */
  /* print_long_array(dependencies, dlength); */

  while(now < duration) {
    /* printf("step %i\n", step); */

    for (up = 0; up < update_length; up++) {
      reaction = update[up];
      propensities[reaction] = rates[reaction];

      for (reactant = 0; reactant < reactants_lengths[reaction]; reactant++) {
        index = reactants_indexes[reaction] + reactant;
        count = outcome[reactants[index]];
        /* printf("count of %ld = %f\n", reactants[index], count); */
        propensities[reaction] *= choose(count, reactions[index]);
      }

      /* printf("updated propensity[%ld]: %f\n", reaction, propensities[reaction]); */
    }

    /* printf("propensities: "); */
    /* print_array(propensities, reactions_length); */

    total = 0.0;
    for (reaction = 0; reaction < reactions_length; reaction++) {
      total += propensities[reaction];
    }

    /* printf("total: %f\n", total); */

    if (total == 0.0) {
      interval = 0.0;
      choice = -1;
      /* printf("breaking because total is 0.0\n"); */
      break;
    } else {
      sample = (double) rand() / RAND_MAX;
      interval = -log(1 - sample) / total;
      point = ((double) rand() / RAND_MAX) * total;

      /* printf("sample: %f - interval: %f - point: %f\n", sample, interval, point); */

      choice = 0;
      progress = 0.0;
      while(progress + propensities[choice] < point) {
        progress += propensities[choice];
        choice += 1;
      }

      /* printf("choice: %d\n", choice); */
      /* printf("progress: %f - now: %f - duration: %f\n", progress, now, duration); */

      if (choice == -1 || (now + interval) > duration) {
        /* printf("completing: %d %d\n", choice == -1, (now + interval) > duration); */
        break;
      }

      now += interval;

      /* printf("new time: %f\n", now); */

      time[step] = now;
      events[step] = choice;

      actors_length = actors_lengths[choice];
      for (actor = 0; actor < actors_length; actor++) {
        index = actors_indexes[choice] + actor;
        adjustment = stoichiometry[choice * substrates_length + actors[index]];
        outcome[actors[index]] += adjustment;

        /* printf("actor[%ld]: %ld -> %ld = %f\n", index, actor, actors[index], adjustment); */
      }

      /* for (species = 0; species < substrates_length; species++) { */
      /*   state[species] += stoichiometry[choice * substrates_length + species]; */
      /* } */

      update_length = dependencies_lengths[choice];
      for (up = 0; up < update_length; up++) {
        index = dependencies_indexes[choice] + up;
        update[up] = dependencies[index];
      }

      /* printf("update\n"); */
      /* print_long_array(update, update_length); */

      step += 1;
    }
  }

  /* print_array(state, substrates_length); */
  
  evolve_result result = {step, time, events, outcome};

  free(propensities);
  free(update);

  return result;
}
