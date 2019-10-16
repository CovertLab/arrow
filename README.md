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

Add the following to your `requirements.txt`, or
`pip install stochastic-arrow`:

    stochastic-arrow

## Usage

The `arrow` library presents a single class as an interface,
`StochasticSystem`, which operates on a set of reactions (encoded as a `numpy`
matrix of stoichiometrix coefficients) and associated reaction rates:

```python
from arrow import StochasticSystem
import numpy as np

# Each row is a reaction and each column is a molecular species (or other
# entity). The first reaction here means that the first and second elements
# combine to create the third, while the fourth is unaffected.
stoichiometric_matrix = np.array([
    [1, 1, -1, 0],
    [-2, 0, 0, 1],
    [-1, -1, 1, 0]], np.int64)

# Each reaction has an associated rate for how probable that reaction is.
rates = np.array([3.0, 1.0, 1.0])

# Once we have a matrix of reactions and their associated rates, we can
# construct the system.
system = StochasticSystem(stoichiometric_matrix, rates)
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
```

Once we have an initial state and duration, we can run the simulation for the
given duration. `evolve` returns a dictionary with five keys:

* steps - the number of steps the simulation took
* time - at what time point each event took place
* events - the events that occurred
* occurrences - the number of times each event occurred (derived directly from `events`)
* outcome - the final state of the system

```python
result = system.evolve(state, duration)
```

If you are interested in the history of states for plotting or otherwise, these can be
derived from the list of events and the stoichiometric matrix, along with the inital
state. `reenact_events` will do this for you:

```python
from arrow import reenact_events

history = reenact_events(stoichiometry, result['events'], state)
```

## Testing

`arrow` uses [pytest](https://docs.pytest.org/en/latest/). To test it:

    > make clean compile
    > pytest

**NOTE:** `make compile` without an explicit `clean` might not fully build the extension.

There are more command line features in test_arrow:

    > python -m arrow.test.test_arrow --complexation

    > python -m arrow.test.test_arrow --plot

    > python -m arrow.test.test_arrow --obsidian

    > python -m arrow.test.test_arrow --memory

    > python -m arrow.test.test_arrow --time
