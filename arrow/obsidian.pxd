# cython: language_level=3str

from libc.stdint cimport int64_t

from mersenne cimport MTState


cdef extern from "obsidian.h":

    ctypedef struct evolve_result:
        int status        # 0 => success
        int steps         # -1 => failure
        double *time      # double time[steps]
        int64_t *events   # int64_t events[steps]
        int64_t *outcome  # int64_t outcome[substrates_count]

    ctypedef struct Info:
        MTState *random_state

        int reactions_count
        int substrates_count
        int64_t *stoichiometry
        # double *rates

        int64_t *reactants_lengths
        int64_t *reactants_indexes
        int64_t *reactants
        int64_t *reactions

        int64_t *dependencies_lengths
        int64_t *dependencies_indexes
        int64_t *dependencies

        int64_t *substrates_lengths
        int64_t *substrates_indexes
        int64_t *substrates

    evolve_result evolve(Info *info, double duration, int64_t *state, double *rates)

    int print_array(double *array, int length)
