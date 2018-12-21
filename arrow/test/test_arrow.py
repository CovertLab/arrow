'''

This submodule is for the expression of small to moderate test systems.  Pytest
will catch any function called 'test_*'.  If called directly, this submodule
will run all systems (provided they are enumerated in the 'systems' module-
level attribute defined below) and plot their output, assuming the correct
return  pattern is used.

'''

from __future__ import absolute_import, division, print_function

import os

import numpy as np
import json

from arrow import evolve, StochasticSystem

def test_equilibration():
    stoichiometric_matrix = np.array([
        [-1, +1,  0],
        [+1, -1, -1],
        ])

    rates = np.array([10, 10, 0.1])
    system = StochasticSystem(stoichiometric_matrix, rates)

    state = np.array([1000, 0])
    duration = 10

    time, counts, events = system.evolve(state, duration)

    assert counts[-1].sum() < state.sum()
    assert time[-1] <= duration

    return (time, counts, events)


def test_dimerization():
    stoichiometric_matrix = np.array([
        [-1, -2, +1],
        [-1,  0, +1],
        [+1,  0, -1],
        [ 0, +1,  0]
        ])

    rates = np.array([3, 1, 1]) * 0.01
    system = StochasticSystem(stoichiometric_matrix, rates)

    state = np.array([1000, 1000, 0, 0])
    duration = 1

    time, counts, events = system.evolve(state, duration)

    assert time[-1] <= duration

    return (time, counts, events)


def test_complexation():
    fixtures_root = os.path.join('data', 'complexation')

    def load_state(filename):
        with open(os.path.join(fixtures_root, filename)) as f:
            state = np.array(json.load(f))

        return state

    initial_state = load_state('initial_state.json')
    final_state = load_state('final_state.json')

    assert initial_state.size == final_state.size

    n_metabolites = initial_state.size

    with open(os.path.join(fixtures_root, 'stoichiometry.json')) as f:
        stoichiometry_sparse = json.load(f)

    n_reactions = len(stoichiometry_sparse)

    stoichiometric_matrix = np.zeros((n_metabolites, n_reactions), np.int64)

    for (reaction_index, reaction_stoich) in enumerate(stoichiometry_sparse):
        for (str_metabolite_index, stoich) in reaction_stoich.viewitems():
            # JSON doesn't allow for integer keys...
            metabolite_index = int(str_metabolite_index)
            stoichiometric_matrix[metabolite_index, reaction_index] = stoich

    duration = 1

    # semi-quantitative rate constants
    rates = np.full(n_reactions, 10)

    system = StochasticSystem(stoichiometric_matrix, rates)

    time, counts, events = system.evolve(initial_state, duration)

    assert(len(time)-1 == events.sum())

    outcome = counts[-1]
    difference = (final_state - outcome)

    total = np.abs(difference).sum()

    print('differences: {}'.format(total))
    print('total steps: {}'.format(len(time)))
    print(time)

    return (time, counts, events)

if __name__ == '__main__':
    from itertools import izip

    import matplotlib.pyplot as plt

    from arrow.analysis.plotting import plot_full_history

    systems = (
        test_equilibration,
        test_dimerization,
        test_complexation,
        )

    n_systems = len(systems)

    ncols = int(np.ceil(np.sqrt(n_systems)))
    nrows = int(np.ceil(n_systems / ncols))

    margins = 1
    axes_size = 3

    figsize = (
        margins + axes_size*ncols,
        margins + axes_size*nrows
        )

    (fig, all_axes) = plt.subplots(
        figsize = figsize,
        nrows = nrows, ncols = ncols,
        constrained_layout = True
        )

    all_axes = np.asarray(all_axes)

    for (axes, system) in izip(all_axes.flatten(), systems):
        axes.set_title(system.func_name)

        time, counts, events = system()
        plot_full_history(axes, time, counts)

    fig.savefig('test_systems.png')
