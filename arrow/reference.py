from __future__ import absolute_import, division, print_function

import numpy as np
from six import moves

from arrow.math import multichoose


def derive_reactants(stoichiometric_matrix):
    '''
    Calculate the various derived values this Gillespie implementation uses extensively
    to avoid recalculating these values every step or call to `evolve`.

    Args:
        stoichiometric_matrix: A matrix representing all of the reactions the system
            is capable of. Each row is a reaction, and each column is a substrate.

    Returns:
        reactants: Array of indexes into each reactant (substrate consumed by the
            reaction) for each reaction.
        reactant_stoichiometries: The value of the reaction for each reactant involved.
        substrates: Array of indexes for each substrate involved in the reaction
            (reactant or product).
    '''

    reactants = [
        np.where(reaction < 0)[0]
        for reaction in stoichiometric_matrix]

    reactant_stoichiometries = [
        -stoichiometric_matrix[reaction, r]
        for reaction, r in enumerate(reactants)]

    substrates = [
        np.where(reaction != 0)[0]
        for reaction in stoichiometric_matrix]

    return reactants, reactant_stoichiometries, substrates

def calculate_dependencies(stoichiometric_matrix):
    '''
    Find out which reactions depend on which other reactions. A dependency exists if
    one of the reactants or products of a reaction is a substrate involved in another reaction.

    Args:
        stoichiometric_matrix: A matrix representing all of the reactions the system
            is capable of. Each row is a reaction, and each column is a substrate.

    Returns:
        dependencies: Array of indexes for each reaction that points to other reactions
            that will be affected by the reaction.
    '''

    dependencies = [
        np.where(np.any(stoichiometric_matrix.T[reaction != 0] < 0, 0))[0]
        for reaction in stoichiometric_matrix]

    return dependencies

def step(
        stoichiometric_matrix,
        rates,
        state,
        reactants=None,
        reactant_stoichiometries=None,
        propensities=None,
        update_reactions=None):
    '''
    Determines which reaction happens next and how much time passes before this event
    given the stoichiometric_matrix, the rates of each reaction, and counts of all
    substrates.
    '''

    if reactants is None:
        reactants, reactant_stoichiometries, substrates = derive_reactants(
            stoichiometric_matrix)

    if update_reactions is None:
        n_reactions = stoichiometric_matrix.shape[0]

        propensities = np.empty(n_reactions)
        update_reactions = moves.range(n_reactions)

    for reaction in update_reactions:
        propensities[reaction] = rates[reaction] * multichoose(
            state[reactants[reaction]],
            reactant_stoichiometries[reaction])

    total = propensities.sum()
    choice = None

    if total == 0:
        interval = 0
        outcome = state

    else:
        interval = np.random.exponential(1 / total)
        random = np.random.uniform(0, 1) * total

        progress = 0
        for choice, span in enumerate(propensities):
            progress += span
            if random <= progress:
                break

        reaction = stoichiometric_matrix[choice]
        outcome = state + reaction

    return interval, outcome, choice, propensities


def evolve(
        stoichiometric_matrix,
        rates,
        state,
        duration,
        reactants=None,
        reactant_stoichiometries=None,
        dependencies=None):
    '''
    Perform a series of steps in the Gillepsie algorithm until the given duration is
    reached.
    '''

    now = 0
    time = [0]
    events = []
    counts = [state]
    propensities = None
    update_reactions = None
    occurrences = np.zeros(rates.shape)

    if reactants is None or reactant_stoichiometries is None:
        reactants, reactant_stoichiometries, substrates = derive_reactants(
            stoichiometric_matrix)

    if dependencies is None:
        dependencies = calculate_dependencies(stoichiometric_matrix)

    while True:
        interval, state, choice, propensities = step(
            stoichiometric_matrix,
            rates,
            state,
            reactants,
            reactant_stoichiometries,
            propensities,
            update_reactions)

        now += interval

        if choice is None or now > duration:
            break

        time.append(now)
        events.append(choice)
        counts.append(state)

        occurrences[choice] += 1

        update_reactions = dependencies[choice]

    time = np.array(time)
    counts = np.array(counts)

    return {
        'steps': len(events),
        'time': np.array(time),
        'events': np.array(events),
        'occurrences': occurrences,
        'outcome': counts[-1],
        'counts': np.array(counts)}

class GillespieReference(object):
    '''
    This is a pure python implementation of the Gillespie algorithm:
    https://en.wikipedia.org/wiki/Gillespie_algorithm

    It has been subsumed by the native implementation in C, Obsidian, but is still
    useful here as reference to the algorithm.
    '''

    def __init__(self, stoichiometric_matrix):
        self.stoichiometric_matrix = stoichiometric_matrix

        reactants, reactant_stoichiometries, substrates = derive_reactants(
            stoichiometric_matrix)

        self.reactants = reactants
        self.reactant_stoichiometries = reactant_stoichiometries
        self.dependencies = calculate_dependencies(stoichiometric_matrix)

    def step(self, state, rates):
        return step(self.stoichiometric_matrix, rates, state)

    def evolve(self, duration, state, rates):
        return evolve(
            self.stoichiometric_matrix,
            rates,
            state,
            duration,
            reactants=self.reactants,
            reactant_stoichiometries=self.reactant_stoichiometries,
            dependencies=self.dependencies)
