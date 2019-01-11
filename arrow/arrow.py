
from __future__ import absolute_import, division, print_function

from itertools import izip

import numpy as np

from arrow.math import multichoose

def _get_reactants_and_stoichiometries(stoichiometric_matrix):
    reactants = [
        np.where(stoichiometry < 0)[0]
        for stoichiometry in stoichiometric_matrix.T
        ]
    reactant_stoichiometries = [
        -stoichiometric_matrix[r, reaction_index]
        for reaction_index, r in enumerate(reactants)
        ]

    return reactants, reactant_stoichiometries


def step(
        stoichiometric_matrix,
        rates,
        state,
        reactants = None, reactant_stoichiometries = None,
        propensities=None, update_reactions=None
        ):

    if reactants is None:
        reactants, reactant_stoichiometries = _get_reactants_and_stoichiometries(
            stoichiometric_matrix
            )

    if update_reactions is None:
        n_reactions = stoichiometric_matrix.shape[1]

        propensities = np.empty(n_reactions)
        update_reactions = xrange(n_reactions)

    for reaction_index in update_reactions:
        propensities[reaction_index] = rates[reaction_index] * multichoose(
            state[reactants[reaction_index]],
            reactant_stoichiometries[reaction_index]
            )

    total = propensities.sum()

    if total == 0:
        time_to_next = 0
        outcome = state
        choice = None

    else:
        time_to_next = np.random.exponential(1 / total)
        random = np.random.uniform(0, 1) * total

        progress = 0
        for choice, interval in enumerate(propensities):
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
    propensities = None
    update_reactions = None
    events = np.zeros(rates.shape)

    reactants, reactant_stoichiometries = _get_reactants_and_stoichiometries(
        stoichiometric_matrix
        )

    dependencies = [
        np.where(np.any(stoichiometric_matrix[stoichiometry != 0] < 0, 0))[0]
        for stoichiometry in stoichiometric_matrix.T]

    while True:
        time_to_next, state, choice, propensities = step(
            stoichiometric_matrix,
            rates,
            state,
            reactants, reactant_stoichiometries,
            propensities, update_reactions
            )

        time_current += time_to_next

        if choice is None or time_current > duration:
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
