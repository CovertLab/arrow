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

    history, steps, events = system.evolve(state, duration)

    assert history[-1].sum() < state.sum()
    assert steps[-1] <= duration

    return (history, steps, events)


def test_dimers():
    reactions = np.array([
        [-1, -1, 1, 0],
        [-2, 0, 0, 1],
        [1, 1, -1, 0]])

    rates = np.array([3, 1, 1]) * 0.01
    system = StochasticSystem(reactions, rates)

    state = np.array([1000, 1000, 0, 0])
    duration = 1

    history, steps, events = system.evolve(state, duration)

    assert steps[-1] <= duration

    return (history, steps, events)


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

    system = StochasticSystem(stoichiometry, rates)

    history, steps, events = system.evolve(state, duration)

    assert(len(steps)-1 == events.sum())

    outcome = history[-1]
    difference = (expected - outcome)
    total = np.abs(difference).sum()

    print('differences: {}'.format(total))
    print('total steps: {}'.format(len(steps)))
    print(steps)

    return (history, steps, events)
    

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
        history, steps, events = system()
        plot_full_history(axes, steps, history)

    plt.savefig('test_systems.png')
