#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "mersenne.h"
#include "obsidian.h"

// Initial length of event vectors
static const int INITIAL_LENGTH = 4000;

// Create a value to represent failures
static const evolve_result failure = {-1, NULL, NULL, NULL};

// Find the number of combinations of choosing k selections from n items
double
choose(long n, long k) {
  double combinations = 1.0;
  long i;
  for (i = 0; i < k; i++) {
    combinations *= ((double) (n - i)) / (i + 1);
  }
  return combinations;
}

// Perform the Gillespie algorithm with the given stoichiometry, reaction rates and
// initial state for the provided duration.

// The main arguments to this function are:
//   * reactions_count: the number of reactions in this system.
//   * substrates_count: the number of substrates this system operates on.
//   * stoichiometry: an array of length `reactions_count * substrates_count` that
//       contains all the information about the reactions this system performs. Each
//       row is a reaction, and each column is a substrate, so every entry is the
//       change in counts of each substrate when the given reaction is performed.
//   * rates: an array of length `reactions_count` that encodes the base rate for
//       each reaction. The actual propensity is further dependent on the counts for
//       each reactant in the reaction.
//   * duration: How long to run the simulation for.
//   * state: An array of length `substrates_count` that contains the inital count
//       for each substrate. The outcome of the algorithm will involve an array of
//       the same size that signifies the counts of each substrate after all the
//       reactions are performed.

// There are four arrays that are derived from the original stoichiometry that are
// nonetheless passed into this function to avoid recomputing these every step. These
// values are nested arrays whose subarrays are of variable length. Instead of messing
// around with pointers to pointers (which would still require an accompanying `lengths`
// array for each nested array), this implementation renders these nested arrays as a
// single flat array with two corresponding arrays describing the index into each
// subarray and the length of each subarray (*_indexes and *_lengths respectively).

// So to iterate through a subarray of array `array` at index `sub`:

//   for (long i = 0; i < lengths[sub]; i++) {
//     long index = indexes[sub];
//     long value = array[index + i];
//     ..... (do something with value)
//   }

// You will see this pattern throughout this function.

// The derived values are:
//   * reactants: Array of indexes into each reactant (substrate consumed by the
//       reaction) for each reaction.
//   * reactions: The value of the reaction for each reactant involved.
//   * substrates: Array of indexes for each substrate involved in the reaction
//       (reactant or product).
//   * dependencies: Array of indexes for each reaction that points to other reactions
//       that will be affected by the reaction.

// The return value of this function is a struct of type `evolve_result`, defined in
// obsidian.h. This has four fields:
//   * steps: How many steps were performed
//   * time: An array of length `steps` containing the value at each time point
//   * events: An array of length `steps` signifying which reaction took place at
//       each time point
//   * outcome: The final state after all of the reactions have been performed.
evolve_result
evolve(MTState *random_state,

       int reactions_count,
       int substrates_count,
       long *stoichiometry,
       double *rates,

       long *reactants_lengths,
       long *reactants_indexes,
       long *reactants,
       long *reactions,
       
       long *dependencies_lengths,
       long *dependencies_indexes,
       long *dependencies,

       long *substrates_lengths,
       long *substrates_indexes,
       long *substrates,

       double duration,
       long *state) {

  // The `event_bounds` will be used to determine how much space to allocate for
  // tracking the evolution of the system's state. If a step is reached that exceeds
  // `event_bounds` then these arrays will be reallocated after doubling `event_bounds`.
  long event_bounds = INITIAL_LENGTH;

  // Allocate the dynamic arrays that will be used to track the progress of the system
  double *time = malloc((sizeof (double)) * event_bounds);
  long *events = malloc((sizeof (long)) * event_bounds);
  long *outcome = malloc((sizeof (long)) * substrates_count);

  // Allocate space for the temporary values that will be used entirely within this
  // function
  double *propensities = malloc((sizeof (double)) * reactions_count);
  long *update = malloc((sizeof (long)) * reactions_count);
  long update_length = reactions_count;

  if (time == NULL ||
      events == NULL ||
      outcome == NULL ||
      propensities == NULL ||
      update == NULL) {
    printf("arrow.obsidian.evolve - failed to allocate memory: %d", errno);

    free(time);
    free(events);
    free(outcome);
    free(propensities);
    free(update);

    return failure;
  }

  // Declare the working variables we will use throughout this function
  long substrates_length;
  long reaction, reactant, index, involve, count, adjustment;
  double total, interval, point, progress;
  int choice, step = 0, up = 0;
  double now = 0.0;

  // Copy the initial state that was supplied from outside to the working `outcome`
  // array we will use to actually apply the reactions and determine the next step's
  // state.
  memcpy(outcome, state, (sizeof (long)) * substrates_count);

  // The `update` array will hold for each step which propensities need to be updated
  // for the next time step, based on the dependencies between reactions (sharing
  // substrates).
  for (reaction = 0; reaction < reactions_count; reaction++) {
    update[reaction] = reaction;
  }

  // Calculate steps until we reach the provided duration
  while (now < duration) {

    // First update all propensities that were affected by the previous event
    for (up = 0; up < update_length; up++) {

      // Find which reaction to update and initialize the propensity for that reaction
      // with the rate for that reaction
      reaction = update[up];
      propensities[reaction] = rates[reaction];

      // Go through each reactant and calculate its contribution to the propensity
      // based on the counts of the corresponding substrate, which is then multiplied
      // by the reaction's original rate and the contributions from other reactants
      for (reactant = 0; reactant < reactants_lengths[reaction]; reactant++) {
        index = reactants_indexes[reaction] + reactant;
        count = outcome[reactants[index]];
        propensities[reaction] *= choose(count, reactions[index]);
      }
    }

    // Find the total for all propensities
    total = 0.0;
    for (reaction = 0; reaction < reactions_count; reaction++) {
      total += propensities[reaction];
    }

    // If the total is zero, then we have no more reactions to perform and can exit
    // early
    if (total == 0.0) {
      interval = 0.0;
      choice = -1;
      break;

    // Otherwise we need to find the next reaction to perform. 
    } else {

      // First, sample two random values, `point` from a linear distribution and
      // `interval` from an exponential distribution.

      interval = sample_exponential(random_state, total);
      point = sample_uniform(random_state);

      // Based on the random sample, find the event that it corresponds to by
      // iterating through the propensities until we surpass our sampled value
      choice = 0;
      progress = 0.0;
      while (progress + propensities[choice] < point) {
        progress += propensities[choice];
        choice += 1;
      }

      // If we have surpassed the provided duration we can exit now
      if (choice == -1 || (now + interval) > duration) {
        break;
      }

      // Increase time by the interval sampled above
      now += interval;

      // Record the information about the chosen event this step
      time[step] = now;
      events[step] = choice;

      // For each substrate involved in this reaction, update the ongoing state
      // with its value from the stoichiometric matrix
      substrates_length = substrates_lengths[choice];
      for (involve = 0; involve < substrates_length; involve++) {
        index = substrates_indexes[choice] + involve;
        adjustment = stoichiometry[choice * substrates_count + substrates[index]];
        outcome[substrates[index]] += adjustment;
      }

      // Find which propensities depend on this reaction and therefore need to be
      // updated in the next step (the rest of the propensities will not change)
      update_length = dependencies_lengths[choice];
      for (up = 0; up < update_length; up++) {
        index = dependencies_indexes[choice] + up;
        update[up] = dependencies[index];
      }

      // Advance
      step += 1;

      // If our step has advanced beyond the current `event_bounds`, double the
      // `event_bounds` and reallocate these arrays with the new size
      if (step >= event_bounds) {
        double *new_time = malloc((sizeof (double)) * event_bounds * 2);
        if (new_time == NULL) {
          printf("arrow.obsidian.evolve - failed to allocate memory: %d", errno);

          free(time);
          free(events);
          free(outcome);
          free(propensities);
          free(update);

          return failure;
        }

        memcpy(new_time, time, (sizeof (double)) * event_bounds);
        free(time);
        time = new_time;

        long *new_events = malloc((sizeof (long)) * event_bounds * 2);
        if (new_events == NULL) {
          printf("arrow.obsidian.evolve - failed to allocate memory: %d", errno);

          free(time);
          free(events);
          free(outcome);
          free(propensities);
          free(update);

          return failure;
        }

        memcpy(new_events, events, (sizeof (long)) * event_bounds);
        free(events);
        events = new_events;

        event_bounds *= 2;
      }
    }
  }

  // Construct the `evolve_result` from the results of performing steps until either
  // the desired duration was achieved or there are no more reactions to perform.
  // We return:
  //   * How many steps were performed
  //   * An array of the time of each event
  //   * An array of each event that occurred each time step
  //   * The resulting state of the system after the chosen reactions were applied
  evolve_result result = {step, time, events, outcome};

  // Clean up the transient allocations
  free(propensities);
  free(update);

  // Return the result constructed above
  return result;
}

// Print an array of doubles
int
print_array(double *array, int length) {
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

// Print an array of longs
int
print_long_array(long *array, int length) {
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
