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


johns_system()
test_dimers()
