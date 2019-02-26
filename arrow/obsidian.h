#ifndef OBSIDIAN_H
#define OBSIDIAN_H

#include "mersenne.h"

// The structure for holding the result of the Gillespie algorithm
typedef struct evolve_result evolve_result;
struct evolve_result {
  int steps;
  double *time;
  int64_t *events;
  int64_t *outcome;
};

// Invoke the system with all the required information to run for the given duration
evolve_result
evolve(MTState *random_state,

       int reactions_count,
       int substrates_count,
       int64_t *stoichiometry,
       double *rates,

       int64_t *reactants_lengths,
       int64_t *reactants_indexes,
       int64_t *reactants,
       int64_t *reactions,
       
       int64_t *dependencies_lengths,
       int64_t *dependencies_indexes,
       int64_t *dependencies,

       int64_t *substrates_lengths,
       int64_t *substrates_indexes,
       int64_t *substrates,

       double duration,
       int64_t *state);

// Supporting print utilities
int print_array(double *array, int length);
int print_int64_t_array(int64_t *array, int length);

#endif
