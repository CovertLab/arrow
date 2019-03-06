#ifndef PARALLEL_H
#define PARALLEL_H

#include "mersenne.h"

typedef struct reaction_state reaction_state;
struct reaction_state {
  int64_t *stoichiometry; // read only
  int64_t total_reactions;
  int64_t total_substrates;
  int64_t index;
  double rate;

  int64_t substrates_length;
  int64_t *substrates;
}

int
evolve_reaction(MTState *random_state,
                reaction_state *reaction,
                double time_delta,
                int64_t *state);

#endif
