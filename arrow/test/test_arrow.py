import numpy as np

from arrow import StochasticSystem


def johns_system():
	reactions = np.array([
		[-1, 1],
		[1, -1],
		[0, -1]])

	weights = np.array([10, 10, 0.1])
	system = StochasticSystem(reactions, weights)

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

	weights = np.array([3, 1, 1])
	system = StochasticSystem(reactions, weights)

	state = np.array([1000, 1000, 0, 0])
	duration = 1
	history, steps = system.evolve(state, duration)

	assert steps[-1] <= duration

	return (history, steps)


johns_system()
test_dimers()
