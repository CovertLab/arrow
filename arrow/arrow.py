import numpy as np


def choose(n, k):
    terms = [
        (n + 1.0 - i) / i
        for i in xrange(1, k + 1)]
    product = np.array(terms).prod()
    return np.rint(product)


def propensity(stoichiometry, state):
    reactants = np.where(stoichiometry < 0)
    terms = [
        choose(state[reactant], -stoichiometry[reactant])
        for reactant in reactants[0]]

    return np.array(terms).prod()


def step(stoichiometric_matrix, rates, state, propensities=[], update_reactions=()):
    if len(update_reactions):
        for update in update_reactions:
            stoichiometry = stoichiometric_matrix[:, update]
            propensities[update] = propensity(stoichiometry, state)
    else:
        propensities = np.array([
            propensity(stoichiometry, state)
            for index, stoichiometry in enumerate(stoichiometric_matrix.T)])

    distribution = (rates * propensities)
    total = distribution.sum()

    if total == 0:
        time_to_next = 0
        outcome = state
        choice = -1

    else:
        time_to_next = np.random.exponential(1 / total)
        random = np.random.uniform(0, 1) * total

        progress = 0
        for choice, interval in enumerate(distribution):
            progress += interval
            if random <= progress:
                break

        stoichiometry = stoichiometric_matrix[:, choice]
        outcome = state + stoichiometry

    return time_to_next, outcome, choice, propensities


def evolve(stoichiometric_matrix, rates, state, duration):
    time_current = 0
    time = [0]
    counts = [state]
    propensities = []
    update_reactions = []
    events = np.zeros(rates.shape)

    dependencies = [
        np.where(np.any(stoichiometric_matrix[stoichiometry != 0] < 0, 0))[0]
        for stoichiometry in stoichiometric_matrix.T]

    while True:
        time_to_next, state, choice, propensities = step(
            stoichiometric_matrix,
            rates,
            state,
            propensities,
            update_reactions)

        time_current += time_to_next
        if not time_to_next or time_current > duration:
            break

        time.append(time_current)
        counts.append(state)

        events[choice] += 1

        update_reactions = dependencies[choice]

    time = np.array(time)
    counts = np.array(counts)

    return time, counts, events


class StochasticSystem(object):
    def __init__(self, stoichiometric_matrix, rates):
        self.stoichiometric_matrix = stoichiometric_matrix
        self.rates = rates

    def step(self, state):
        return step(self.stoichiometric_matrix, self.rates, state)

    def evolve(self, state, duration):
        return evolve(self.stoichiometric_matrix, self.rates, state, duration)
