import numpy as np


def choose(n, k):
    terms = [
        (n + 1.0 - i) / i
        for i in xrange(1, k + 1)]
    product = np.array(terms).prod()
    return np.rint(product)


def propensity(stoichiometric_matrix, state, form):
    reactants = np.where(stoichiometric_matrix < 0)
    terms = [
        form(state[reactant], -stoichiometric_matrix[reactant])
        for reactant in reactants[0]]

    return np.array(terms).prod()

def step(stoichiometric_matrix, rates, state, forms, propensities=[], update_reactions=()):
    if len(update_reactions):
        for update in update_reactions:
            stoichiometry = stoichiometric_matrix[:, update]
            form = forms if callable(forms) else forms[update]
            propensities[update] = propensity(stoichiometry, state, form)
    else:
        propensities = np.array([
            propensity(stoichiometry, state, forms if callable(forms) else forms[index])
            for index, stoichiometry in enumerate(stoichiometric_matrix.T)])

    distribution = (rates * propensities)
    total = distribution.sum()

    if total == 0:
        dt = 0
        outcome = state
        choice = -1

    else:
        dt = np.random.exponential(1 / total)
        random = np.random.uniform(0, 1) * total

        progress = 0
        for choice, interval in enumerate(distribution):
            progress += interval
            if random <= progress:
                break

        stoichiometry = stoichiometric_matrix[:, choice]
        outcome = state + stoichiometry

    return dt, outcome, choice, propensities


def evolve(stoichiometric_matrix, rates, state, duration, forms=choose):
    t = 0
    time = [0]
    counts = [state]
    propensities = []
    update_reactions = []

    while True:
        dt, state, choice, propensities = step(
            stoichiometric_matrix,
            rates,
            state,
            forms,
            propensities,
            update_reactions)

        t += dt
        if not dt or t > duration:
            break

        counts.append(state)
        time.append(t)

        stoichiometry = stoichiometric_matrix[:, choice]
        involved = np.where(stoichiometry != 0)[0]
        update_reactions = np.where(np.any(stoichiometric_matrix[involved] < 0, 0))[0]

    time = np.array(time)
    counts = np.array(counts)

    return time, counts


class StochasticSystem(object):
    def __init__(self, stoichiometric_matrix, rates, forms=None):
        self.stoichiometric_matrix = stoichiometric_matrix
        self.rates = rates
        self.forms = forms or choose

    def step(self, state):
        return step(self.stoichiometric_matrix, self.rates, state, forms=self.forms)

    def evolve(self, state, duration):
        return evolve(self.stoichiometric_matrix, self.rates, state, duration, forms=self.forms)
