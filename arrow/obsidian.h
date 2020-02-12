#ifndef OBSIDIAN_H
#define OBSIDIAN_H

#include "mersenne.h"

// The structure for holding the result of the Gillespie algorithm
typedef struct evolve_result {
  int steps;         // -1 => failure
  double *time;      // double time[steps]
  int64_t *events;   // int64_t events[steps]
  int64_t *outcome;  // int64_t outcome[substrates_count]
} evolve_result;

typedef struct Info {
    MTState *random_state;

    int reactions_count;
    int substrates_count;
    int64_t *stoichiometry;

    int64_t *reactants_lengths;
    int64_t *reactants_indexes;
    int64_t *reactants;
    int64_t *reactions;

    int64_t *dependencies_lengths;
    int64_t *dependencies_indexes;
    int64_t *dependencies;

    int64_t *substrates_lengths;
    int64_t *substrates_indexes;
    int64_t *substrates;
} Info;

// Invoke the system with all the required information to run for the given duration.
// The result is either a failure = {-1, NULL, NULL, NULL} or it points to malloc'd
// arrays that the caller must free().
evolve_result evolve(Info *info, double duration, int64_t *state, double *rates);

// Supporting print utilities
int print_array(double *array, int length);

#endif
