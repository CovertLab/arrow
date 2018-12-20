'''

This submodule is for the expression of small to moderate test systems.  Pytest
will catch any function called 'test_*'.  If called directly, this submodule
will run all systems (provided they are enumerated in the 'systems' module-
level attribute defined below) and plot their output, assuming the correct
return  pattern is used.

'''

from __future__ import absolute_import, division, print_function

import numpy as np
import json

from arrow import evolve, StochasticSystem

def test_equilibration():
    reactions = np.array([
        [-1, +1,  0],
        [+1, -1,  0],
        [ 0,  0, -1]
        ])

    rates = np.array([10, 10, 0.1])
    system = StochasticSystem(reactions, rates)

    state = np.array([1000, 0])
    duration = 1

    time, counts = system.evolve(state, duration)

    assert counts[-1].sum() < state.sum()
    assert time[-1] <= duration

    return (time, counts)


def test_dimerization():
    reactions = np.array([
        [-1, -2, +1],
        [-1,  0, +1],
        [+1,  0, -1],
        [ 0, +1,  0]
        ])

    rates = np.array([3, 1, 1]) * 0.01
    system = StochasticSystem(reactions, rates)

    state = np.array([1000, 1000, 0, 0])
    duration = 1

    time, counts = system.evolve(state, duration)

    assert time[-1] <= duration

    return (time, counts)


def test_complexation():
    with open('data/complexation/complexation-state.json', 'r') as file:
        data = json.load(file)

    assert len(data) == 3

    stoichiometry = np.array(data['stoichiometry']).transpose()
    state = np.array(data['before'])
    expected = np.array(data['after'])
    duration = 1

    assert len(state) == len(expected)
    assert stoichiometry.shape[1] == len(state)

    # semi-quantitative rate constants
    rates = np.full(stoichiometry.shape[0], 10)

    system = StochasticSystem(stoichiometry.T, rates)

    time, counts = system.evolve(state, duration)

    outcome = counts[-1]
    difference = (expected - outcome)
    total = np.abs(difference).sum()

    print('differences: {}'.format(total))
    print('total steps: {}'.format(len(time)))
    print(time)

    return (time, counts)


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

    for (axes, system) in izip(all_axes.flatten(), systems):
        axes.set_title(system.func_name)
        plot_full_history(axes, *system())

    fig.savefig('test_systems.png')
