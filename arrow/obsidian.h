int print_array(double * array, int length);

typedef struct evolve_result evolve_result;
struct evolve_result {
  int steps;
  double * time;
  int * events;
  double * state;
};

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
       double * state);
