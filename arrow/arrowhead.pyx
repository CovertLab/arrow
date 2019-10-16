# cython: language_level=3str

from __future__ import absolute_import, division, print_function

cimport cython
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from libc.stdint cimport int64_t
from libc.string cimport memset, memcpy
from libc.stdlib cimport free

import numpy as np
cimport numpy as np

from . cimport mersenne
from . cimport obsidian


np.import_array()  # Initialize numpy's C API so it won't segfault


cdef class Arrowhead:
    """Cython interface to the C-coded Gillespie algorithm "obsidian", leaving
    obsidian.c free of Python dependencies.
    """

    cdef:
        int random_seed
        obsidian.Info info
        tuple refs

    # __cinit__() is always called immediately on construction, before CPython
    # considers calling __init__(), which it can skip. `self` is not yet fully
    # constructed so *BE CAREFUL* invoking Python operations that might touch it!
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __cinit__(self, *args):
        memset(&self.info, 0, sizeof(obsidian.Info))
        self.info.random_state = <mersenne.MTState*> PyMem_Malloc(sizeof(mersenne.MTState))
        if not self.info.random_state:
            raise MemoryError()

    # This will accept any type of array in C-contiguous layout that support
    # the memoryview protocol (numpy array, array.array, Cython array, buffer
    # interface) with the declared element type and number of dimensions.
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __init__(
            self,
            int random_seed,
            int64_t[:, ::1] stoichiometry not None,
            double[::1] rates not None,

            int64_t[::1] reactants_lengths not None,
            int64_t[::1] reactants_indexes not None,
            int64_t[::1] reactants not None,
            int64_t[::1] reactions not None,

            int64_t[::1] dependencies_lengths not None,
            int64_t[::1] dependencies_indexes not None,
            int64_t[::1] dependencies not None,

            int64_t[::1] substrates_lengths not None,
            int64_t[::1] substrates_indexes not None,
            int64_t[::1] substrates not None,
            *args):
        self.random_seed = random_seed
        mersenne.seed(self.info.random_state, random_seed)

        self.info.reactions_count = stoichiometry.shape[0]
        self.info.substrates_count = stoichiometry.shape[1]
        self.info.stoichiometry = &stoichiometry[0, 0]
        self.info.rates = &rates[0]

        self.info.reactants_lengths = &reactants_lengths[0]
        self.info.reactants_indexes = &reactants_indexes[0]
        self.info.reactants = &reactants[0]
        self.info.reactions = &reactions[0]

        self.info.dependencies_lengths = &dependencies_lengths[0]
        self.info.dependencies_indexes = &dependencies_indexes[0]
        self.info.dependencies = &dependencies[0]

        self.info.substrates_lengths = &substrates_lengths[0]
        self.info.substrates_indexes = &substrates_indexes[0]
        self.info.substrates = &substrates[0]

        # TODO(jerry): Check array sizes? E.g. rates.shape[0] >= reactions_count

        self.refs = (  # hold refs to these memoryviews while we hold their data ptrs
            stoichiometry, rates,
            reactants_lengths, reactants_indexes, reactants, reactions,
            dependencies_lengths, dependencies_indexes, dependencies,
            substrates_lengths, substrates_indexes, substrates)

    # `self` might no longer be a valid Python object so *BE CAREFUL* invoking
    # Python operations that might touch it!
    def __dealloc__(self):
        PyMem_Free(self.info.random_state)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def evolve(self, double duration, int64_t[::1] state):
        """Run the Gillespie algorithm with the initialized stoichiometry and
        rates over the given duration and state.
        Return None or a tuple (steps, time, events, outcome).
        """
        # TODO(jerry): Check the state[] array size?

        evolved = obsidian.evolve(&self.info, duration, &state[0])
        cdef int steps = evolved.steps
        cdef int count = self.info.substrates_count

        if steps == -1:
            return None

        time = copy_c_array(evolved.time, steps, sizeof(double), np.NPY_DOUBLE)
        events = copy_c_array(evolved.events, steps, sizeof(int64_t), np.NPY_INT64)
        outcome = copy_c_array(evolved.outcome, count, sizeof(int64_t), np.NPY_INT64)

        result = steps, time, events, outcome

        free(evolved.time)
        free(evolved.events)
        free(evolved.outcome)

        return result

    def reactions_count(self):
        """Returns the number of reactions for this system."""
        return self.info.reactions_count

    def substrates_count(self):
        """Returns the number of substrates this system operates on."""
        return self.info.substrates_count

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def demo(self):
        """A demo that prints the reaction rates to C stdout."""
        obsidian.print_array(self.info.rates, self.info.reactions_count)

cdef np.ndarray copy_c_array(
        void *source, np.npy_intp element_count, size_t element_size, int np_typenum):
    """Copy the source 1-D C array into a numpy array so it's independent of the
    source buffer.
    """
    cdef np.ndarray result = np.PyArray_SimpleNew(1, &element_count, np_typenum)
    memcpy(np.PyArray_DATA(result), source, element_size * element_count)
    return result
