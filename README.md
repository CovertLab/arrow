# Arrow

“... even if the previous millisecond is closer to us than the birth of the universe, it is equally out of reach.”
― Jean-Christophe Valtat, Luminous Chaos

## Concept

This library implements a generalized version of the [Gillespie
Algorithm](https://en.wikipedia.org/wiki/Gillespie_algorithm), a stochastic
approach to numerically solving discrete systems. Each iteration, the algorithm
will calculate the propensities for each reaction given a rate and the counts
of the reactants present in the current state of the system, then selects one
reaction to occur and the interval of time between the previous reaction and
the current reaction. Iterating this produces a trajectory (or `history`) of
the state vector over the course of the simulation.

## Installation

Add the following to your `requirements.txt`, or `pip install stochastic-
arrow`:

    stochastic-arrow==0.0.1

## Usage

The `arrow` library presents a single class as an interface,
`StochasticSystem`, which operates on a set of reactions (encoded as a `numpy`
matrix of stoichiometrix coefficients) and associated reaction rates:

```python
from arrow import StochasticSystem
import numpy as np

# Each column is a reaction and each row is a molecular species (or other
# entity). The first reaction here means that the first and second elements
# combine to create the third, while the fourth is unaffected.
stoichiometry = np.array([
    [-1, -2, +1],
    [-1,  0, +1],
    [+1,  0, -1],
    [ 0, +1,  0]
    ])

# Each reaction has an associated rate for how probable that reaction is.
rates = np.array([3.0, 1.0, 1.0])

# Once we have a matrix of reactions and their associated rates, we can
# construct the system.
system = StochasticSystem(stoichiometry, rates)
```

Now that the system has been instantiated, we can invoke it with any initial
state vector and then run it for a given time interval:

```python
# This gives the initial state of the system (counts of each molecular species,
# for instance).
state = np.array([1000, 1000, 0, 0])

# We also specify how long we want the simulation to run. Here we set it to one
# second.
duration = 1

# Once we have an initial state and duration, we can run the simulation for the
# given duration. `evolve` returns the history of the state vector (counts) for
# each time step, and the history of time steps as they will be in uneven
# increments throughout the simulation.
time, counts = system.evolve(state, duration)
```

## Testing

`arrow` uses pytest: https://docs.pytest.org/en/latest/ so you can test simply
by invoking:

    > pytest

Also, we have a test that generates plots of various systems which can be run
like so:

    > python arrow/test/test_arrow.py
