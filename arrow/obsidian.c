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

static int INITIAL_LENGTH = 4000;

double
choose(long n, long k) {
  double combinations = 1.0;
  long i;
  for (i = 0; i < k; i++) {
    combinations *= ((double) (n - i)) / (i + 1);
  }
  return combinations;
}

evolve_result
evolve(int reactions_length,
       int substrates_length,
       long * stoichiometry,
       double * rates,

       long * reactants_lengths,
       long * reactants_indexes,
       long * reactants,
       long * reactions,
       
       long * dependencies_lengths,
       long * dependencies_indexes,
       long * dependencies,

       long * actors_lengths,
       long * actors_indexes,
       long * actors,

       double duration,
       long * state) {

  long event_bounds = INITIAL_LENGTH;
  double * time = malloc((sizeof (double *)) * event_bounds);
  long * events = malloc((sizeof (long *)) * event_bounds);
  long * outcome = malloc((sizeof (long *)) * substrates_length);

  double * propensities = malloc((sizeof (double *)) * reactions_length);
  long * update = malloc((sizeof (long *)) * reactions_length);
  long update_length = reactions_length;
  long actors_length;

  long reaction, reactant, species, index, actor, count, adjustment;
  double total, interval, sample, progress, point;
  
  int choice, step = 0, up = 0;
  double now = 0.0;

  for (species = 0; species < substrates_length; species++) {
    outcome[species] = state[species];
  }

  for (reaction = 0; reaction < reactions_length; reaction++) {
    update[reaction] = reaction;
  }

  while(now < duration) {
    for (up = 0; up < update_length; up++) {
      reaction = update[up];
      propensities[reaction] = rates[reaction];

      for (reactant = 0; reactant < reactants_lengths[reaction]; reactant++) {
        index = reactants_indexes[reaction] + reactant;
        count = outcome[reactants[index]];
        propensities[reaction] *= choose(count, reactions[index]);
      }
    }

    total = 0.0;
    for (reaction = 0; reaction < reactions_length; reaction++) {
      total += propensities[reaction];
    }

    if (total == 0.0) {
      interval = 0.0;
      choice = -1;
      break;
    } else {
      sample = (double) rand() / RAND_MAX;
      interval = -log(1 - sample) / total;
      point = ((double) rand() / RAND_MAX) * total;

      choice = 0;
      progress = 0.0;
      while(progress + propensities[choice] < point) {
        progress += propensities[choice];
        choice += 1;
      }

      if (choice == -1 || (now + interval) > duration) {
        break;
      }

      now += interval;

      time[step] = now;
      events[step] = choice;

      actors_length = actors_lengths[choice];
      for (actor = 0; actor < actors_length; actor++) {
        index = actors_indexes[choice] + actor;
        adjustment = stoichiometry[choice * substrates_length + actors[index]];
        outcome[actors[index]] += adjustment;
      }

      update_length = dependencies_lengths[choice];
      for (up = 0; up < update_length; up++) {
        index = dependencies_indexes[choice] + up;
        update[up] = dependencies[index];
      }

      step += 1;

      if (step >= event_bounds) {
        double * new_time = malloc((sizeof (double *)) * event_bounds * 2);
        memcpy(new_time, time, (sizeof (double *)) * event_bounds);
        free(time);
        time = new_time;

        long * new_events = malloc((sizeof (long *)) * event_bounds * 2);
        memcpy(new_events, events, (sizeof (long *)) * event_bounds);
        free(events);
        events = new_events;

        event_bounds *= 2;
      }
    }
  }

  evolve_result result = {step, time, events, outcome};

  free(propensities);
  free(update);

  return result;
}
