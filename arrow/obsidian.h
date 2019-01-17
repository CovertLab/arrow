int print_array(double * array, int length);
int print_long_array(long * array, int length);

typedef struct evolve_result evolve_result;
struct evolve_result {
  int steps;
  double * time;
  long * events;
  long * outcome;
};

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
       long * state);
