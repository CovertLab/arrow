# cython: language_level=3str
# cython: freethreading_compatible = True

from libc.stdint cimport uint32_t


cdef extern from "mersenne.h":

    size_t TWISTER_SIZE

    ctypedef struct MTState:
        # mersenne.h defines TWISTER_SIZE to be 624.
        uint32_t MT[624]
        uint32_t MT_TEMPERED[624]
        size_t index


    void seed(MTState *state, uint32_t seed_value)

    uint32_t rand_u32(MTState *state)

    double sample_uniform(MTState *state)

    double sample_exponential(MTState *state, double lambda_)
