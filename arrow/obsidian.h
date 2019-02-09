#ifndef OBSIDIAN_H
#define OBSIDIAN_H

#include "mersenne.h"

// The structure for holding the result of the Gillespie algorithm
typedef struct evolve_result evolve_result;
struct evolve_result {
  int steps;
  double *time;
  long *events;
  long *outcome;
};

// Invoke the system with all the required information to run for the given duration
evolve_result
evolve(MTState *random_state,

       int reactions_count,
       int substrates_count,
       long *stoichiometry,
       double *rates_flat,
       long *rates_lengths,
       long *rates_indexes,
       long *forms,

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
       long *state);

// Supporting print utilities
int print_array(double *array, int length);
int print_long_array(long *array, int length);

#endif
