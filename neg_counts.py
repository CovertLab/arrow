import json

from arrow import StochasticSystem
from arrow.reference import GillespieReference
import numpy as np


duration = 2**31

with open('neg_counts.json') as f:
    data = json.load(f)

stoich = np.array(data['stoich'])
rates = np.array(data['rates'])
counts = np.array(data['counts'])

rates = np.full(rates.shape, 1.0)

system = GillespieReference(stoich.T, rates)
# system = StochasticSystem(stoich.T, rates, random_seed=0)

while True:
    result = system.evolve(duration, counts)

    updated_counts = result['outcome']

    if not np.any(counts - updated_counts):
        break

    if np.any(updated_counts < 0):
        negative = np.where(updated_counts < 0)[0]
        print('negative indexes: {}'.format(negative))
        print('negative counts: {}'.format(updated_counts[negative]))

        raise Exception('Negative counts')
    
    counts = updated_counts
    print(counts)
 
