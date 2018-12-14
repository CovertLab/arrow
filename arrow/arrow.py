import numpy as np
from scipy.special import comb
from numba import jit


@jit
def choose(pair):
	reaction = pair[0]
	molecules = pair[1]
	return comb(molecules, -reaction)

@jit
def propensity(reaction, state):
	reactants = np.where(reaction < 0)
	terms = map(choose, zip(reaction[reactants], state[reactants]))
	return np.array(terms).prod()

@jit
def step(reactions, rates, state, propensities=[], update_reactions=()):
	if len(update_reactions):
		for update in update_reactions:
			reaction = reactions[update]
			propensities[update] = propensity(reaction, state)
	else:
		propensities = [
			propensity(reaction, state)
			for reaction in reactions]

	distribution = (rates * propensities)
	total = distribution.sum()
	if total == 0:
		return (state, 0, -1, propensities)

	dt = np.random.exponential(1 / total)
	random = np.random.uniform(0, 1) * total
		
	progress = 0
	for choice, interval in enumerate(distribution):
		progress += interval
		if random <= progress:
			break

	reaction = reactions[choice]
	outcome = state + reaction

	return outcome, dt, choice, propensities

@jit
def evolve(reactions, rates, state, duration):
	time = 0
	history = [state]
	steps = [0]
	propensities = []
	update_reactions = []

	while True:
		state, dt, choice, propensities = step(
			reactions,
			rates,
			state,
			propensities,
			update_reactions)

		time += dt
		if time > duration:
			break

		history.append(state)
		steps.append(time)

		reaction = reactions[choice]
		involved = np.where(reaction != 0)
		update_reactions = np.where(reactions[:, involved] != 0)[0]

	return history, steps


class StochasticSystem(object):
	def __init__(self, reactions, rates):
		self.reactions = reactions
		self.rates = rates

	def step(self, state):
		return step(self.reactions, self.rates, state)

	def evolve(self, state, duration):
		return evolve(self.reactions, self.rates, state, duration)
