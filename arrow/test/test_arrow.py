'''

This submodule is for the expression of small to moderate test systems.  Pytest
will catch any function called 'test_*'.  If called directly, this submodule
will run all systems (provided they are enumerated in the 'systems' module-
level attribute defined below) and plot their output, assuming the correct
return  pattern is used.

'''

from __future__ import absolute_import, division, print_function

from six import moves
import os
from time import time as seconds_since_epoch
import json
import numpy as np
import psutil
import argparse
import pickle

from arrow import reenact_events, StochasticSystem
from arrow import GillespieReference

def test_equilibration():
    stoichiometric_matrix = np.array([
        [0, -1],
        [+1, -1],
        [-1, +1]])

    rates = np.array([10, 10, 0.1])
    system = GillespieReference(stoichiometric_matrix)

    state = np.array([1000, 0])
    duration = 10

    result = system.evolve(duration, state, rates)

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
    system = GillespieReference(stoichiometric_matrix)

    state = np.array([1000, 1000, 0, 0])
    duration = 1

    result = system.evolve(duration, state, rates)

    time = result['time']
    counts = reenact_events(stoichiometric_matrix, result['events'], state)
    events = result['occurrences']

    assert time[-1] <= duration

    return (time, counts, events)


def load_complexation(prefix='simple'):
    fixtures_root = os.path.join('data', 'complexation')

    def load_state(filename):
        with open(os.path.join(fixtures_root, filename)) as f:
            state = np.array(json.load(f))

        return state

    initial_state = load_state(prefix + '-initial.json')
    final_state = load_state(prefix + '-final.json')

    assert initial_state.size == final_state.size

    n_metabolites = initial_state.size

    with open(os.path.join(fixtures_root, 'stoichiometry.json')) as f2:
        stoichiometry_sparse = json.load(f2)

    n_reactions = len(stoichiometry_sparse)

    stoichiometric_matrix = np.zeros((n_reactions, n_metabolites), np.int64)

    for (reaction_index, reaction_stoich) in enumerate(stoichiometry_sparse):
        for (str_metabolite_index, stoich) in reaction_stoich.items():
            # JSON doesn't allow for integer keys...
            metabolite_index = int(str_metabolite_index)
            stoichiometric_matrix[reaction_index, metabolite_index] = stoich

    # semi-quantitative rate constants
    rates = np.full(n_reactions, 1000.0)

    return (stoichiometric_matrix, rates, initial_state, final_state)

def complexation_test(make_system):
    stoichiometric_matrix, rates, initial_state, final_state = load_complexation(prefix='simple')
    duration = 1

    system = make_system(stoichiometric_matrix)
    result = system.evolve(duration, initial_state, rates)

    time = np.concatenate([[0.0], result['time']])
    events = result['events']
    occurrences = result['occurrences']
    outcome = result['outcome']

    history = reenact_events(stoichiometric_matrix, events, initial_state)

    difference = (final_state - outcome)

    total = np.abs(difference).sum()

    print('counts before: {}'.format(initial_state.sum()))
    print('counts after: {}'.format(outcome.sum()))
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

    arrow = StochasticSystem(stoichiometric_matrix)
    result = arrow.evolve(1.0, np.array([50, 20, 30, 40], np.int64), rates)

    print('steps: {}'.format(result['steps']))
    print('time: {}'.format(result['time']))
    print('events: {}'.format(result['events']))
    print('occurrences: {}'.format(result['occurrences']))
    print('outcome: {}'.format(result['outcome']))

    assert(arrow.obsidian.reactions_count() == stoichiometric_matrix.shape[0])
    assert(arrow.obsidian.substrates_count() == stoichiometric_matrix.shape[1])

    return result

def test_compare_runtime():
    stoichiometric_matrix, rates, initial_state, final_state = load_complexation()
    duration = 1
    amplify = 100

    reference = GillespieReference(stoichiometric_matrix)
    reference_start = seconds_since_epoch()
    for i in range(amplify):
        result = reference.evolve(duration, initial_state, rates)
    reference_end = seconds_since_epoch()

    system = StochasticSystem(stoichiometric_matrix)
    obsidian_start = seconds_since_epoch()
    for i in range(amplify):
        result = system.evolve(duration, initial_state, rates)
    obsidian_end = seconds_since_epoch()

    print('reference Python implementation elapsed seconds: {}'.format(
        reference_end - reference_start))
    print('obsidian C implementation elapsed seconds: {}'.format(
        obsidian_end - obsidian_start))

def test_memory():
    stoichiometric_matrix, rates, initial_state, final_state = load_complexation()
    duration = 1
    amplify = 100

    this = psutil.Process(os.getpid())
    memory = memory_previous = this.memory_info().rss
    memory_increases = 0
    print('initial memory use: {}'.format(memory))

    system = StochasticSystem(stoichiometric_matrix, random_seed=np.random.randint(2**31))

    obsidian_start = seconds_since_epoch()
    for i in range(1, amplify + 1):
        memory = this.memory_info().rss
        if memory > memory_previous:
            print('memory use before iteration {:2d}: {}'.format(i, memory))
            memory_previous = memory
            memory_increases += 1

        result = system.evolve(duration, initial_state, rates)
        difference = np.abs(final_state - result['outcome']).sum()

        if difference:
            print('difference is {}'.format(difference))
    obsidian_end = seconds_since_epoch()

    print('obsidian C implementation elapsed seconds for {} runs: {}'.format(
        amplify, obsidian_end - obsidian_start))
    assert memory_increases <= 1

def test_pickle():
    stoichiometric_matrix = np.array([
        [1, 1, -1, 0],
        [-2, 0, 0, 1],
        [-1, -1, 1, 0]], np.int64)

    rates = np.array([3, 1, 1]) * 0.01

    arrow = StochasticSystem(stoichiometric_matrix)

    pickled_arrow = pickle.dumps(arrow)
    unpickled_arrow = pickle.loads(pickled_arrow)

    result = arrow.evolve(1.0, np.array([50, 20, 30, 40], np.int64), rates)

    straight = test_obsidian()

    assert(result['steps'] == straight['steps'])
    assert((result['time'] == straight['time']).all())
    assert((result['events'] == straight['events']).all())
    assert((result['occurrences'] == straight['occurrences']).all())
    assert((result['outcome'] == straight['outcome']).all())

    print('arrow object pickled is {} bytes'.format(len(pickled_arrow)))


def main(args):
    systems = (
        test_equilibration,
        test_dimerization,
        lambda: complexation_test(StochasticSystem),
        lambda: complexation_test(GillespieReference))

    if not args.plot:
        if args.complexation:
            for run in range(args.runs):
                complexation_test(StochasticSystem)
        if args.obsidian:
            test_obsidian()
        elif args.memory:
            test_memory()
        elif args.time:
            test_compare_runtime()
        elif args.pickle:
            test_pickle()
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
            margins + axes_size * ncols,
            margins + axes_size * nrows)

        (fig, all_axes) = plt.subplots(
            figsize=figsize,
            nrows=nrows, ncols=ncols,
            constrained_layout=True)

        all_axes = np.asarray(all_axes)

        for (axes, system) in moves.zip(all_axes.flatten(), systems):
            axes.set_title(system.func_name)

            time, counts, events = system()
            plot_full_history(axes, time, counts)

        fig.savefig('test_systems.png')
        print('Wrote test_systems.png')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Run one of these tests')
    parser.add_argument('--complexation', action='store_true')
    parser.add_argument('--runs', type=int, default=1)
    parser.add_argument('--plot', action='store_true')
    parser.add_argument('--obsidian', action='store_true')
    parser.add_argument('--memory', action='store_true')
    parser.add_argument('--time', action='store_true')
    parser.add_argument('--pickle', action='store_true')

    main(parser.parse_args())
