'''

This submodule is for the expression of small to moderate test systems.  Pytest
will catch any function called 'test_*'.  If called directly, this submodule
will run all systems (provided they are enumerated in the 'systems' module-
level attribute defined below) and plot their output, assuming the correct
return  pattern is used.

'''

from __future__ import absolute_import, division, print_function

import os
import time
import json
import numpy as np
import argparse

from arrow import reenact_events, StochasticSystem
from arrow import derive_reactants, calculate_dependencies, evolve, GillespieReference

import obsidian

def test_equilibration():
    stoichiometric_matrix = np.array([
        [0, -1],
        [+1, -1],
        [-1, +1]])

    rates = np.array([10, 10, 0.1])
    system = GillespieReference(stoichiometric_matrix, rates)

    state = np.array([1000, 0])
    duration = 10

    result = system.evolve(duration, state)

    time = result['time']
    counts = reenact_events(stoichiometric_matrix, result['events'], state)
    events = result['occurrences']

    assert counts[-1].sum() < state.sum()
    assert time[-1] <= duration

    return (time, counts, events)


def test_dimerization():
    stoichiometric_matrix = np.array([
        [1, 1, -1, 0],
        [-2, 0, 0, 1],
        [-1, -1, 1, 0]], np.int64)

    rates = np.array([3, 1, 1]) * 0.01
    system = GillespieReference(stoichiometric_matrix, rates)

    state = np.array([1000, 1000, 0, 0])
    duration = 1

    result = system.evolve(duration, state)

    time = result['time']
    counts = reenact_events(stoichiometric_matrix, result['events'], state)
    events = result['occurrences']

    assert time[-1] <= duration

    return (time, counts, events)


def load_complexation():
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

    stoichiometric_matrix = np.zeros((n_reactions, n_metabolites), np.int64)

    for (reaction_index, reaction_stoich) in enumerate(stoichiometry_sparse):
        for (str_metabolite_index, stoich) in reaction_stoich.viewitems():
            # JSON doesn't allow for integer keys...
            metabolite_index = int(str_metabolite_index)
            stoichiometric_matrix[reaction_index, metabolite_index] = stoich

    duration = 1

    # semi-quantitative rate constants
    rates = np.full(n_reactions, 1000)

    return stoichiometric_matrix, rates, initial_state, final_state

def complexation_test(make_system):
    stoichiometric_matrix, rates, initial_state, final_state = load_complexation()
    duration = 1

    system = make_system(stoichiometric_matrix, rates)
    result = system.evolve(duration, initial_state)

    time = np.concatenate([[0.0], result['time']])
    events = result['events']
    occurrences = result['occurrences']
    outcome = result['outcome']

    history = reenact_events(stoichiometric_matrix, events, initial_state)

    difference = (final_state - outcome)

    total = np.abs(difference).sum()

    print('differences: {}'.format(total))
    print('total steps: {}'.format(len(time)))
    print('number of events: {}'.format(len(events)))
    print('length of history: {}'.format(len(history)))
    print(time)

    return (time, history, occurrences)

def test_complexation():
    complexation_test(StochasticSystem)

def test_obsidian():
    stoichiometric_matrix = np.array([
        [1, 1, -1, 0],
        [-2, 0, 0, 1],
        [-1, -1, 1, 0]], np.int64)

    rates = np.array([3, 1, 1]) * 0.01

    arrow = StochasticSystem(stoichiometric_matrix, rates)
    result = arrow.evolve(1.0, np.array([50, 20, 30, 40], np.int64))

    print('steps: {}'.format(result['steps']))
    print('time: {}'.format(result['time']))
    print('events: {}'.format(result['events']))
    print('occurrences: {}'.format(result['occurrences']))
    print('outcome: {}'.format(result['outcome']))

    assert(arrow.obsidian.reactions_count() == stoichiometric_matrix.shape[0])
    assert(arrow.obsidian.substrates_count() == stoichiometric_matrix.shape[1])

def test_compare_runtime():
    stoichiometric_matrix, rates, initial_state, final_state = load_complexation()
    duration = 1
    amplify = 100

    reference = GillespieReference(stoichiometric_matrix, rates)
    reference_start = time.time()
    for i in range(amplify):
        result = reference.evolve(duration, initial_state)
    reference_end = time.time()

    system = StochasticSystem(stoichiometric_matrix, rates)
    obsidian_start = time.time()
    for i in range(amplify):
        result = system.evolve(duration, initial_state)
    obsidian_end = time.time()

    print('reference time elapsed: {}'.format(reference_end - reference_start))
    print('obsidian time elapsed: {}'.format(obsidian_end - obsidian_start))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--plot', action='store_true')
    parser.add_argument('--complexation', action='store_true')
    parser.add_argument('--runs', type=int, default=1)
    args = parser.parse_args()

    from itertools import izip

    systems = (
        test_equilibration,
        test_dimerization,
        lambda: complexation_test(StochasticSystem),
        lambda: complexation_test(GillespieReference))

    if not args.plot:
        if args.complexation:
            for run in xrange(args.runs):
                lambda: complexation_test(StochasticSystem)
        else:
            for system in systems:
                system()
    else:
        import matplotlib.pyplot as plt
        from arrow.analysis.plotting import plot_full_history

        n_systems = len(systems)

        ncols = int(np.ceil(np.sqrt(n_systems)))
        nrows = int(np.ceil(n_systems / ncols))

        margins = 1
        axes_size = 3

        figsize = (
            margins + axes_size*ncols,
            margins + axes_size*nrows)

        (fig, all_axes) = plt.subplots(
            figsize = figsize,
            nrows = nrows, ncols = ncols,
            constrained_layout = True)

        all_axes = np.asarray(all_axes)

        for (axes, system) in izip(all_axes.flatten(), systems):
            axes.set_title(system.func_name)

            time, counts, events = system()
            plot_full_history(axes, time, counts)

        fig.savefig('test_systems.png')
