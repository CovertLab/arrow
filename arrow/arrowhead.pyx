# cython: language_level=3str

from cpython.mem cimport PyMem_Malloc, PyMem_Free

import numpy as np
cimport numpy as np

from . cimport mersenne
from . cimport obsidian


cdef class Arrowhead:
    """Cython interface to the C-coded Gillespie algorithm."""

    cdef:
        int random_seed
        obsidian.Info info

    # __cinit__ is always called immediately on construction, before CPython
    # considers calling __init__ which can be skipped. self is not fully
    # constructed so tread lightly on it besides initializing cdef fields.
    def __cinit__(self, int random_seed):
        self.info.random_state = self.info.stoichiometry = self.info.rates = \
            self.info.reactants_lengths = self.info.reactants_indexes = \
            self.info.reactants = \
            self.info.reactions = self.info.dependencies_lengths = \
            self.info.dependencies_indexes = self.info.dependencies = \
            self.info.substrates_lengths = self.info.substrates_indexes = \
            self.info.substrates = NULL

        self.random_seed = random_seed
        self.info.random_state = <mersenne.MTState*> PyMem_Malloc(sizeof(mersenne.MTState))
        if not self.info.random_state:
            raise MemoryError()
        mersenne.seed(self.info.random_state, random_seed)

    def __dealloc__(self):
        PyMem_Free(self.info.random_state)
