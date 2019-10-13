# cython: language_level=3str

from libc.stdint cimport int64_t

from mersenne cimport MTState


cdef extern from "obsidian.h":

    ctypedef struct evolve_result:
        pass

    evolve_result evolve(
        MTState *random_state,

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
        int64_t *state)

    int print_array(double *array, int length)
