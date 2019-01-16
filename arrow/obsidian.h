int print_array(double * array, int length);

typedef struct evolve_result evolve_result;
struct evolve_result {
  double * time;
  double ** counts;
  int * events;
};

evolve_result
evolve(int reactions_length,
       int substrates_length,
       double * stoichiometry,
       double * rates,
       double * state,
       double duration);
