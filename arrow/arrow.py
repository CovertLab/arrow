from __future__ import absolute_import, division, print_function

import numpy as np
from .reference import derive_reactants, calculate_dependencies
from .arrowhead import Arrowhead


def flatten(l):
    '''
    Flatten a list by one level:
        [[1, 2, 3], [[4, 5], 6], [7]] --> [1, 2, 3, [4, 5], 6, 7]
    '''

    return [
        item
        for sublist in l
        for item in sublist]

def flat_indexes(assorted_lists):
    '''
    Take a list of variable length lists and reduce it to a single flat array with
    associated indexes and lengths. These indexes and lengths can be used in
    conjunction with the flat array to recover the original list of lists.

    Args:
        assorted_lists (List[List]): A list of variable length lists. 

    Returns numpy arrays:
        flat: The flattened data.
        lengths: The lengths of all the given sublists.
        indexes: The indexes where the sublists begin in the flat array.
    '''

    lengths = np.array([
        len(l)
        for l in assorted_lists])
    indexes = np.insert(lengths, 0, 0).cumsum()[:-1]
    flat = np.array(flatten(assorted_lists))
    return flat, lengths, indexes

def reenact_events(stoichiometry, events, state):
    '''
    Reproduce the history of states given an initial state, the history of events and
    the stoichiometry of those events.

    Args:
        stoichiometry: 2d array where each row is a reaction and each column is a substrate.
        events: History of reactions as indexes into the stoichiometric matrix.
        state: Initial state of the system.

    The function iteratively applies each reaction from `events` to the original state
    and returns the history of states the system took on.

    Return:
        history: 2d array where each row is the state of the system at some point in
            the history of its evolution. Each subsequent row is the previous state with
            the next reaction from `events` applied.
    '''

    history = np.zeros((events.shape[0] + 1, state.shape[0]))
    history[0] = state
    for index, event in enumerate(events):
        history[index + 1] = history[index] + stoichiometry[event]
    return history

class StochasticSystem(object):
    '''
    This class encapsulates a stochastic system which performs the Gillespie algorithm
    with a given stoichiometric matrix representing the reactions of the system:

        https://en.wikipedia.org/wiki/Gillespie_algorithm

    A basic summary is that the algorithm is given a stoichiometric matrix representing
    each reaction that can occur in the system, the rates of each reaction, an initial
    state and a duration, and it finds how much time passes and which reaction occurs
    next. By applying this iteratively you can "evolve" the system through time, one
    reaction at a time.

    The stoichiometric matrix has a reaction for each row, with the values in that row
    encoding how many of each substrate are either consumed or produced by the reaction
    (and zero everywhere else). 
    '''

    def __init__(self, stoichiometry, random_seed=0):
        '''
        This invokes the Obsidian C code (via arrowhead.pyx) with the
        stoichiometry, reaction rates and a variety of derived values. Once constructed,
        this can be invoked by calling `evolve` with a duration and initial state, since
        the stoichiometry will be shared among all calls.

        There are four derived values, each of which is a list of variable length
        lists. In order to pass this into C, these nested lists are flattened and two
        associated lists are constructed: one to hold the indexes for each sublist and
        one to hold the lengths of these sublists. So, to access one of the sublists
        at index `i`, you can find the start index of the sublist inside of `_flat`
        with `_indexes[i]`, then iterate through it for the number of elements given
        by `_lengths[i]`. Strictly speaking this could have been implemented with only
        the lengths (or the indexes, they are in a way duals of each other), but that
        would require additional computation. Since this effort is motivated by
        performance, in order to trade time for space the most expedient approach is
        chosen, even if it is more verbose.
        '''

        self.stoichiometry = stoichiometry.copy()
        self.random_seed = random_seed

        reactants, reactions, substrates = derive_reactants(stoichiometry)
        dependencies = calculate_dependencies(stoichiometry)

        self.reactants_flat, self.reactants_lengths, self.reactants_indexes = flat_indexes(reactants)
        self.reactions_flat = np.array(flatten(reactions))
        self.dependencies_flat, self.dependencies_lengths, self.dependencies_indexes = flat_indexes(
            dependencies)
        self.substrates_flat, self.substrates_lengths, self.substrates_indexes = flat_indexes(substrates)

        self.obsidian = Arrowhead(
            self.random_seed,
            self.stoichiometry,
            self.reactants_lengths,
            self.reactants_indexes,
            self.reactants_flat,
            self.reactions_flat,
            self.dependencies_lengths,
            self.dependencies_indexes,
            self.dependencies_flat,
            self.substrates_lengths,
            self.substrates_indexes,
            self.substrates_flat)

    def evolve(self, duration, state, rates):
        steps, time, events, outcome = self.obsidian.evolve(duration, state, rates)
        occurrences = np.zeros(len(rates))
        for event in events:
            occurrences[event] += 1

        return {
            'steps': steps,
            'time': time,
            'events': events,
            'occurrences': occurrences,
            'outcome': outcome}
