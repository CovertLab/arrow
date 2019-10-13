# cython: language_level=3str

from libc.stdint cimport int64_t
from cpython.mem cimport PyMem_Malloc, PyMem_Free

import numpy as np
cimport numpy as np

from . cimport mersenne
# from obsidian cimport evolve


cdef class Arrowhead:
    """Cython interface to the obsidian C code."""

    cdef:
        int random_seed
        mersenne.MTState *random_state

        int reactions_count
        int substrates_count
        int64_t *stoichiometry
        double *rates

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

    # __cinit__ is always called immediately on construction, before CPython
    # considers calling __init__ which can be skipped. self is not fully
    # constructed so tread lightly on it besides initializing cdef fields.
    def __cinit__(self, int random_seed):
        self.random_state = self.stoichiometry = self.rates = \
            self.reactants_lengths = self.reactants_indexes = self.reactants = \
            self.reactions = self.dependencies_lengths = \
            self.dependencies_indexes = self.dependencies = \
            self.substrates_lengths = self.substrates_indexes = \
            self.substrates = NULL

        self.random_seed = random_seed
        self.random_state = <mersenne.MTState*> PyMem_Malloc(sizeof(mersenne.MTState))
        if not self.random_state:
            raise MemoryError()
        mersenne.seed(self.random_state, random_seed)

    def __dealloc__(self):
        PyMem_Free(self.random_state)
