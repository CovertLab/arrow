from __future__ import absolute_import, division, print_function

import numpy as np
from .reference import derive_reactants, calculate_dependencies
import obsidian

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
    '''
    This class encapsulates a stochastic system which performs the Gillespie algorithm with
    a given stoichiometric matrix representing the reactions of the system:

        https://en.wikipedia.org/wiki/Gillespie_algorithm

    The stoichiometric matrix has a reaction for each row, with the values in that row encoding
    how many of each substrate are either consumed or produced by the reaction (and zero everywhere
    else). 
    '''

    def __init__(self, stoichiometry, rates):
		'''
		This invokes the Obsidian C module with the stoichiometry, reaction rates and a variety of
		derived values. Once constructed, this can be invoked by calling `evolve` with a duration
		and initial state, since the stoichiometry will be shared among all calls.

		There are four derived values, each of which is a list of variable length lists. In order
		to pass this into C, these nested lists are flattened and two associated lists are
		constructed: one to hold the indexes for each sublist and one to hold the lengths of these
		sublists. So, to access one of the sublists at index `i`, you can find the start index of
		the sublist inside of `_flat` with `_indexes[i]`, then iterate through it for the number
		of elements given by `_lengths[i]`. Strictly speaking this could have been implemented
		with only the lengths, but that would require iterating through each length to discover
		the start index of the sublist in the flattened form and since this effort is in order to
		trade time for space the most expedient approach is chosen.
		'''

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
