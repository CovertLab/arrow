from __future__ import absolute_import, division, print_function

from itertools import izip

import numpy as np

from arrow.math import multichoose
import obsidian

def derive_reactants(stoichiometry):
    '''
    Calculate the various derived values this Gillespie implementation uses extensively
    to avoid recalculating these values every step or call to `evolve`.

    Args:
        stoichiometry: A matrix representing all of the reactions the system is capable of.
            Each row is a reaction, and each column is a substrate.

    Returns:
        reactants: Array of indexes into each reactant (substrate consumed by the reaction)
            for each reaction.
        reactions: The value of the reaction for each reactant involved.
        involved: Array of indexes for each substrate involved in the reaction (reactant or product).
    '''

    reactants = [
        np.where(reaction < 0)[0]
        for reaction in stoichiometry]

    reactions = [
        -stoichiometry[reaction, r]
        for reaction, r in enumerate(reactants)]

    involved = [
        np.where(reaction != 0)[0]
        for reaction in stoichiometry]

    return reactants, reactions, involved

def calculate_dependencies(stoichiometry):
    '''
    Find out which reactions depend on which other reactions. A dependency exists if one of the
    reactants or products of a reaction is a substrate involved in another reaction.

    Args:
        stoichiometry: A matrix representing all of the reactions the system is capable of.
            Each row is a reaction, and each column is a substrate.

    Returns:
        dependencies: Array of indexes for each reaction that points to other reactions that will
            be affected by the reaction.
    '''

    dependencies = [
        np.where(np.any(stoichiometry.T[reaction != 0] < 0, 0))[0]
        for reaction in stoichiometry]

    return dependencies

def step(
        stoichiometry,
        rates,
        state,
        reactants=None,
        reactions=None,
        propensities=None,
        update_reactions=None):
    '''
    Determines which reaction happens next and how much time passes before this event
    given the stoichiometry, the rates of each reaction, and counts of all substrates.
    '''

    if reactants is None:
        reactants, reactions, involved = derive_reactants(
            stoichiometry)

    if update_reactions is None:
        n_reactions = stoichiometry.shape[0]

        propensities = np.empty(n_reactions)
        update_reactions = xrange(n_reactions)

    for reaction in update_reactions:
        propensities[reaction] = rates[reaction] * multichoose(
            state[reactants[reaction]],
            reactions[reaction])

    total = propensities.sum()

    if total == 0:
        interval = 0
        outcome = state
        choice = None

    else:
        interval = np.random.exponential(1 / total)
        random = np.random.uniform(0, 1) * total

        progress = 0
        for choice, span in enumerate(propensities):
            progress += span
            if random <= progress:
                break

        reaction = stoichiometry[choice]
        outcome = state + reaction

    return interval, outcome, choice, propensities


def evolve(
        stoichiometry,
        rates,
        state,
        duration,
        reactants=None,
        reactions=None,
        dependencies=None):
    '''
    Perform a series of steps in the Gillepsie algorithm until the given duration is reached.
    '''

    now = 0
    time = [0]
    counts = [state]
    propensities = None
    update_reactions = None
    events = np.zeros(rates.shape)

    if reactants is None or reactions is None:
        reactants, reactions, involved = derive_reactants(
            stoichiometry)

    if dependencies is None:
        dependencies = calculate_dependencies(stoichiometry)

    while True:
        interval, state, choice, propensities = step(
            stoichiometry,
            rates,
            state,
            reactants,
            reactions,
            propensities,
            update_reactions)

        now += interval

        if choice is None or now > duration:
            break

        time.append(now)
        counts.append(state)

        events[choice] += 1

        update_reactions = dependencies[choice]

    time = np.array(time)
    counts = np.array(counts)

    return time, counts, events

class GillespieReference(object):
    '''
    This is a pure python implementation of the Gillespie algorithm:
    https://en.wikipedia.org/wiki/Gillespie_algorithm

    It has been subsumed by the native implementation in C, Obsidian, but is still useful here as
    reference to the algorithm.
    '''
    def __init__(self, stoichiometry, rates):
        self.stoichiometry = stoichiometry
        self.rates = rates

        reactants, reactions, involved = derive_reactants(stoichiometry)

        self.reactants = reactants
        self.reactions = reactions
        self.dependencies = calculate_dependencies(stoichiometry)

    def step(self, state):
        return step(self.stoichiometry, self.rates, state)

    def evolve(self, state, duration):
        return evolve(
            self.stoichiometry,
            self.rates,
            state,
            duration,
            reactants=self.reactants,
            reactions=self.reactions,
            dependencies=self.dependencies)

def flatten(l):
    '''
    Flatten a list by one level: [[1, 2, 3], [[4, 5], 6], [7]] --> [1, 2, 3, [4, 5], 6, 7]
    '''
    return [
        item
        for sublist in l
        for item in sublist]

def flat_indexes(variable):
    '''
    Take a list of variable length lists and reduce it to a single flat list with associated
    lengths and indexes, for use by the obsidian C module.
    '''

    lengths = np.array([
        len(l)
        for l in variable])
    indexes = np.insert(lengths, 0, 0).cumsum()[:-1]
    flat = np.array(flatten(variable))
    return flat, lengths, indexes

def reenact_events(stoichiometry, events, state):
    '''
    Reproduce the history of states given an initial state, the history of events and the
    stoichiometry of those events.
    '''
    history = np.zeros((events.shape[0] + 1, state.shape[0]))
    history[0] = state
    for index, event in enumerate(events):
        history[index + 1] = history[index] + stoichiometry[event]
    return history

class StochasticSystem(object):
    def __init__(self, stoichiometry, rates):
        self.stoichiometry = stoichiometry
        self.rates = rates

        reactants, reactions, involved = derive_reactants(stoichiometry)
        dependencies = calculate_dependencies(stoichiometry)

        self.reactants_flat, self.reactants_lengths, self.reactants_indexes = flat_indexes(reactants)
        self.reactions_flat = np.array(flatten(reactions))
        self.dependencies_flat, self.dependencies_lengths, self.dependencies_indexes = flat_indexes(
            dependencies)
        self.involved_flat, self.involved_lengths, self.involved_indexes = flat_indexes(involved)

        self.obsidian = obsidian.obsidian(
            self.stoichiometry,
            self.rates,
            self.reactants_lengths,
            self.reactants_indexes,
            self.reactants_flat,
            self.reactions_flat,
            self.dependencies_lengths,
            self.dependencies_indexes,
            self.dependencies_flat,
            self.involved_lengths,
            self.involved_indexes,
            self.involved_flat)

    def evolve(self, duration, state):
        steps, time, events, outcome = self.obsidian.evolve(duration, state)
        occurrences = np.zeros(len(self.rates))
        for event in events:
            occurrences[event] += 1

        return {
            'steps': steps,
            'time': time,
            'events': events,
            'occurrences': occurrences,
            'outcome': outcome}
