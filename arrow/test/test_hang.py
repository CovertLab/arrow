# A test case for Issue #48.
#
# This code should reproduce the error: the program hangs after it prints "7100"
# and eventually runs out of memory. You'll have to Ctrl+Z to stop it.
# Ctrl+C won't work.
#
# It hangs only with certain seeds and numbers of molecules. The system can
# evolve with the same number of molecule counts for 7179 iterations before it
# hangs. Adding 1 to all of the molecules causes it to hang at an earlier
# iteration.
#
# TODO: Debug this. Is it caused by a Gillespie algorithm blowup (see below),
# integer overflow Undefined Behavior in C, or something else?
#
# The Gillespie algorithm is prone to explode [symptom?] under certain
# conditions if the exponent term in the choice calculation is too large.
#
# The workaround is to find the problematic reaction and decompose the
# stoichiometry into an equivalent problem with more steps.
#
# The flagella example had something like 170 identical subunits which caused
# the problem. Breaking it into 2+ equivalent reactions fixed it.
#
# It'd be good for the Arrow code to catch this problem when/before it happens
# and at least identify which reactions are problematic.

import os

from arrow import StochasticSystem
import numpy as np


def np_load(filename):
    filepath = os.path.join(os.path.dirname(__file__), filename)
    return np.load(filepath)


def test_hang():
    # TODO: Use a pytest plug-in to timeout after some threshold.

    seed = 807952948

    stoich = np_load('stoich.npy')
    mol = np_load('complex-counts.npy')
    rates = np_load('rates.npy')

    system = StochasticSystem(stoich, random_seed=seed)
    for i in range(10000):
        if i % 100 == 0:
            print(i)

        result = system.evolve(1, mol, rates)


if __name__ == '__main__':
    test_hang()
