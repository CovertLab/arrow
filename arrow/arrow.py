from __future__ import absolute_import, division, print_function

from itertools import izip

import numpy as np

from arrow.math import multichoose
import obsidian

def derive_reactants(stoichiometric_matrix):
    reactants = [
        np.where(stoichiometry < 0)[0]
        for stoichiometry in stoichiometric_matrix.T]

    reactant_stoichiometries = [
        -stoichiometric_matrix[r, reaction_index]
        for reaction_index, r in enumerate(reactants)]

    involved = [
        np.where(stoichiometry != 0)[0]
        for stoichiometry in stoichiometric_matrix.T]

    return reactants, reactant_stoichiometries, involved

def calculate_dependencies(stoichiometric_matrix):
    dependencies = [
        np.where(np.any(stoichiometric_matrix[stoichiometry != 0] < 0, 0))[0]
        for stoichiometry in stoichiometric_matrix.T]

    return dependencies

def step(
        stoichiometric_matrix,
        rates,
        state,
        reactants=None,
        reactant_stoichiometries=None,
        propensities=None,
        update_reactions=None):

    if reactants is None:
        reactants, reactant_stoichiometries, involved = derive_reactants(
            stoichiometric_matrix.T)

    if update_reactions is None:
        n_reactions = stoichiometric_matrix.shape[0]

        propensities = np.empty(n_reactions)
        update_reactions = xrange(n_reactions)

    for reaction_index in update_reactions:
        propensities[reaction_index] = rates[reaction_index] * multichoose(
            state[reactants[reaction_index]],
            reactant_stoichiometries[reaction_index])

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

        stoichiometry = stoichiometric_matrix[choice, :]
        outcome = state + stoichiometry

    return interval, outcome, choice, propensities


def evolve(
        stoichiometric_matrix,
        rates,
        state,
        duration,
        reactants=None,
        reactant_stoichiometries=None,
        dependencies=None):
    now = 0
    time = [0]
    counts = [state]
    propensities = None
    update_reactions = None
    events = np.zeros(rates.shape)

    if reactants is None or reactant_stoichiometries is None:
        reactants, reactant_stoichiometries, involved = derive_reactants(
            stoichiometric_matrix.T)

    if dependencies is None:
        dependencies = calculate_dependencies(stoichiometric_matrix.T)

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

        reactants, reactant_stoichiometries, involved = derive_reactants(stoichiometric_matrix.T)

        self.reactants = reactants
        self.reactant_stoichiometries = reactant_stoichiometries
        self.dependencies = calculate_dependencies(stoichiometric_matrix.T)

    def step(self, state):
        return step(self.stoichiometric_matrix, self.rates, state)

    def evolve(self, state, duration):
        return evolve(
            self.stoichiometric_matrix,
            self.rates,
            state,
            duration,
            reactants=self.reactants,
            reactant_stoichiometries=self.reactant_stoichiometries,
            dependencies=self.dependencies)

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

def flatten(l):
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

class Arrow(object):
    def __init__(self, stoichiometry, rates):
        self.stoichiometry = stoichiometry
        self.rates = rates

        reactants, reactions, involved = derive_reactants(stoichiometry.T)
        dependencies = calculate_dependencies(stoichiometry.T)

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
