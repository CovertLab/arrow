import numpy as np

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

	rates = np.array([3, 1, 1])
	system = StochasticSystem(reactions, rates)

	state = np.array([1000, 1000, 0, 0])
	duration = 1

	history, steps = system.evolve(state, duration)

	assert steps[-1] <= duration

	return (history, steps)


if __name__ == '__main__':
	from itertools import izip

	import matplotlib.pyplot as plt

	from arrow.analysis.plotting import plot_full_history

	systems = (johns_system, test_dimers)

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

	for (axes, system) in izip(all_axes, systems):
		axes.set_title(system.func_name)
		plot_full_history(axes, *system()[::-1])

	plt.savefig('test_systems.png')
