
import os

import numpy as np
import json

from arrow import evolve, StochasticSystem

def johns_system():
    reactions = np.array([
        [-1, 1],
        [1, -1],
        [0, -1]])

    rates = np.array([10, 10, 0.1])
    system = StochasticSystem(reactions, rates)

    state = np.array([1000, 0])
    duration = 1

    history, steps = system.evolve(state, duration)

    assert history[-1].sum() < state.sum()
    assert steps[-1] <= duration

    return (history, steps)


def test_dimers():
    reactions = np.array([
        [-1, -1, 1, 0],
        [-2, 0, 0, 1],
        [1, 1, -1, 0]])

    rates = np.array([3, 1, 1]) * 0.01
    system = StochasticSystem(reactions, rates)

    state = np.array([1000, 1000, 0, 0])
    duration = 1

    history, steps = system.evolve(state, duration)

    assert steps[-1] <= duration

    return (history, steps)


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

    stoichiometry = np.zeros((n_metabolites, n_reactions), np.int64)

    for (reaction_index, reaction_stoich) in enumerate(stoichiometry_sparse):
        for (str_metabolite_index, stoich) in reaction_stoich.viewitems():
            # JSON doesn't allow for integer keys...
            metabolite_index = int(str_metabolite_index)
            stoichiometry[metabolite_index, reaction_index] = stoich

    duration = 1

    # semi-quantitative rate constants
    rates = np.full(stoichiometry.shape[1], 10)

    system = StochasticSystem(stoichiometry.T, rates)

    history, steps = system.evolve(initial_state, duration)

    outcome = history[-1]
    difference = (final_state - outcome)
    total = np.abs(difference).sum()

    print('differences: {}'.format(total))
    print('total steps: {}'.format(len(steps)))
    print(steps)

    return (history, steps)


if __name__ == '__main__':
    from itertools import izip

    import matplotlib.pyplot as plt

    from arrow.plotting import plot_full_history

    systems = (johns_system, test_dimers, test_complexation, test_complexation)

    n_systems = len(systems)

    ncols = int(np.ceil(np.sqrt(n_systems)))
    nrows = int(np.ceil(n_systems / ncols))

    margins = 1
    axes_size = 3

    figsize = (
        margins + axes_size*ncols,
        margins + axes_size*nrows
        )

    (fix, all_axes) = plt.subplots(
        figsize = figsize,
        nrows = nrows, ncols = ncols,
        constrained_layout = True
        )

    for (axes, system) in izip(all_axes.flatten(), systems):
        axes.set_title(system.func_name)
        plot_full_history(axes, *system()[::-1])

    plt.savefig('test_systems.png')
