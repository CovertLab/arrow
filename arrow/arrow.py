import numpy as np


def choose(n, k):
    terms = [
        (n + 1.0 - i) / i
        for i in xrange(1, k + 1)]
    product = np.array(terms).prod()
    return np.rint(product)


def propensity(stoichiometry, state, form):
    reactants = np.where(stoichiometry < 0)
    terms = [
        form(state[reactant], -stoichiometry[reactant])
        for reactant in reactants[0]]

    return np.array(terms).prod()


def step(stoichiometry, rates, state, forms, propensities=[], update_reactions=()):
    if len(update_reactions):
        for update in update_reactions:
            reaction_stoichoiometry = stoichiometry[update]
            form = forms if callable(forms) else forms[update]
            propensities[update] = propensity(reaction_stoichoiometry, state, form)
    else:
        propensities = np.array([
            propensity(reaction_stoichoiometry, state, forms if callable(forms) else forms[index])
            for index, reaction_stoichoiometry in enumerate(stoichiometry)])

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

        reaction_stoichoiometry = stoichiometry[choice]
        outcome = state + reaction_stoichoiometry

    return dt, outcome, choice, propensities


def evolve(stoichiometry, rates, state, duration, forms=choose):
    t = 0
    time = [0]
    counts = [state]
    propensities = []
    update_reactions = []

    while True:
        dt, state, choice, propensities = step(
            stoichiometry,
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

        reaction = stoichiometry[choice]
        involved = np.where(reaction != 0)
        update_reactions = np.where(stoichiometry[:, involved] != 0)[0]

    time = np.array(time)
    counts = np.array(counts)

    return time, counts


class StochasticSystem(object):
    def __init__(self, stoichiometry, rates, forms=None):
        self.stoichiometry = stoichiometry
        self.rates = rates
        self.forms = forms or choose

    def step(self, state):
        return step(self.stoichiometry, self.rates, state, forms=self.forms)

    def evolve(self, state, duration):
        return evolve(self.stoichiometry, self.rates, state, duration, forms=self.forms)
