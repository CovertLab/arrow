# cython: language_level=3str

from libc.stdint cimport uint32_t


cdef extern from "mersenne.h":

    ctypedef struct MTState:
        pass

    size_t TWISTER_SIZE

    void seed(MTState *state, uint32_t seed_value)

    uint32_t rand_u32(MTState *state)

    double sample_uniform(MTState *state)

    double sample_exponential(MTState *state, double lambda_)
