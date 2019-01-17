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

    actors = [
        np.where(stoichiometry != 0)[0]
        for stoichiometry in stoichiometric_matrix.T]

    return reactants, reactant_stoichiometries, actors

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
        reactants, reactant_stoichiometries, actors = derive_reactants(
            stoichiometric_matrix)

    if update_reactions is None:
        n_reactions = stoichiometric_matrix.shape[1]

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

        stoichiometry = stoichiometric_matrix[:, choice]
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
        reactants, reactant_stoichiometries, actors = derive_reactants(
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

        reactants, reactant_stoichiometries, actors = derive_reactants(stoichiometric_matrix)

        self.reactants = reactants
        self.reactant_stoichiometries = reactant_stoichiometries
        self.dependencies = calculate_dependencies(stoichiometric_matrix)

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

def flatten(l):
    return [item for sublist in l for item in sublist]

class Arrow(object):
    def __init__(self, stoichiometry, rates):
        self.stoichiometry = stoichiometry
        self.rates = rates

        reactants, reactions, actors = derive_reactants(stoichiometry.T)
        dependencies = calculate_dependencies(stoichiometry.T)

        self.reactants_lengths = np.array([
            len(reactant)
            for reactant in reactants])
        self.reactants_indexes = np.insert(self.reactants_lengths, 0, 0).cumsum()[:-1]
        self.reactants_flat = np.array(flatten(reactants))
        self.reactions_flat = np.array(flatten(reactions))

        self.dependencies_lengths = np.array([
            len(dependency)
            for dependency in dependencies])
        self.dependencies_indexes = np.insert(self.dependencies_lengths, 0, 0).cumsum()[:-1]
        self.dependencies_flat = np.array(flatten(dependencies))

        self.actors_lengths = np.array([
            len(actor)
            for actor in actors])
        self.actors_indexes = np.insert(self.actors_lengths, 0, 0).cumsum()[:-1]
        self.actors_flat = np.array(flatten(actors))

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
            self.actors_lengths,
            self.actors_indexes,
            self.actors_flat)

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
