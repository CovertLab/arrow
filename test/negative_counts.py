import json
import os

from stochastic_arrow import StochasticSystem
import numpy as np


duration = 2**31

# Build the path to the data file relative to the script's directory
data_path = os.path.join(os.path.dirname(__file__), 'data', 'complexation', 'large-initial.json')

with open(data_path) as f:
    data = json.load(f)

stoich = np.array(data['stoich'])
rates = np.array(data['rates']) * 1e-30
counts = np.array(data['counts'])

system = StochasticSystem(stoich.T, random_seed=0)

while True:
    result = system.evolve(duration, counts, rates)

    updated_counts = result['outcome']

    if not np.any(counts - updated_counts):
        break

    if np.any(updated_counts < 0):
        raise Exception('Negative counts')
    
    counts = updated_counts
    print(counts)
