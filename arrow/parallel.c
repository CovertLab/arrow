#include "mersenne.h"
#include "parallel.h"

// Find the number of combinations of choosing k selections from n items
double
choose(int64_t n, int64_t k) {
  double combinations = 1.0;
  int64_t i;
  for (i = 0; i < k; i++) {
    combinations *= ((double) (n - i)) / (i + 1);
  }
  return combinations;
}

int
choose_reaction(MTState *random_state,
                reaction_state *reaction,
                double time_delta,
                int64_t *state) {
  int64_t choice = 0,
    substrate,
    substrate_index,
    stoichiometry;

  double sample = sample_uniform(random_state);

  int64_t reaction_index = reaction->index * reaction->total_reactions;

  // find the propensity given the counts of each reactant
  double propensity = reaction->rate;
  for (substrate = 0; substrate < reaction->substrates_length; substrate++) {
    substrate_index = reaction->substrates[substrate];
    stoichiometry = reaction->stoichiometry[reaction_index + substrate_index];
    if (stoichiometry < 0) { // we have a reactant
      propensity *= choose(state[substrate_index], -stoichiometry)
    }
  }

  
  // choose whether this reaction occurs based on (propensity * time_delta)
  if (sample < propensity * time_delta) {
    choice = 1;
  }

  return choice;
}
