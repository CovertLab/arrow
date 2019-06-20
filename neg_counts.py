import json

from arrow import StochasticSystem
import numpy as np


duration = 2**31

with open('neg_counts.json') as f:
    data = json.load(f)

stoich = np.array(data['stoich'])
rates = np.array(data['rates'])
counts = np.array(data['counts'])

system = StochasticSystem(stoich.T, rates, random_seed=0)

while True:
    result = system.evolve(duration, counts)

    updated_counts = result['outcome']

    if not np.any(counts - updated_counts):
        break

    if np.any(updated_counts < 0):
        raise Exception('Negative counts')
    
    counts = updated_counts
    print(counts)
