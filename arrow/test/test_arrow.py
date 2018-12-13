import numpy as np

from arrow import StochasticSystem


def johns_system():
	reactions = np.array([
		[-1, 1],
		[1, -1],
		[0, -1]])

	weights = np.array([10, 10, 0.1])

	state = np.array([1000, 0])

	system = StochasticSystem(reactions, weights)
	history, steps = system.evolve(state, 1)

	print(zip(steps, history))


def test_dimers():
	reactions = np.array([
		[-1, -1, 1, 0],
		[-2, 0, 0, 1],
		[1, 1, -1, 0]])

	weights = np.array([3, 1, 1])

	state = np.array([1000, 1000, 0, 0])

	system = StochasticSystem(reactions, weights)
	history, steps = system.evolve(state, 1)

	print(zip(steps, history))
