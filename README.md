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

Add the following to your `requirements.txt`, or run
`pip install stochastic-arrow` to install it [from PyPI](https://pypi.org/project/stochastic-arrow/):

    stochastic-arrow

**NOTE:** If upgrading from a version older than 1.0.0, check if the [`arrow`](https://github.com/arrow-py/arrow) datetime package is installed. If so, uninstall `arrow` before upgrading `stochastic-arrow`, then reinstall `arrow`.

    > pip show arrow
    > pip uninstall arrow
    > pip install stochastic-arrow
    > pip install arrow

## Usage

The `stochastic_arrow` library presents a single class as an interface,
`StochasticSystem`, which operates on a set of reactions (encoded as a `numpy`
matrix of stoichiometrix coefficients) and associated reaction rates:

```python
from stochastic_arrow import StochasticSystem
import numpy as np

# Each row is a reaction and each column is a molecular species (or other
# entity). The first reaction here means that the first and second elements
# combine to create the third, while the fourth is unaffected.
stoichiometric_matrix = np.array([
    [1, 1, -1, 0],
    [-2, 0, 0, 1],
    [-1, -1, 1, 0]], np.int64)

# Once we have a matrix of reactions, we can
# construct the system.
system = StochasticSystem(stoichiometric_matrix)
```

Now that the system has been instantiated, we can invoke it with any initial
state vector and set of reaction rates and then run it for a given time interval:

```python
# This gives the initial state of the system (counts of each molecular species,
# for instance).
import numpy as np
state = np.array([1000, 1000, 0, 0])

# We also specify how long we want the simulation to run. Here we set it to one
# second.
duration = 1

# Each reaction has an associated rate for how probable that reaction is.
rates = np.array([3.0, 1.0, 1.0])
```

Once we have an initial state and rates, we can run the simulation for the
given duration. `evolve` returns a dictionary with five keys:

* steps - the number of steps the simulation took
* time - at what time point each event took place
* events - the events that occurred
* occurrences - the number of times each event occurred (derived directly from `events`)
* outcome - the final state of the system

```python
result = system.evolve(duration, state, rates)
```

If you are interested in the history of states for plotting or otherwise, these can be
derived from the list of events and the stoichiometric matrix, along with the inital
state. `reenact_events` will do this for you:

```python
from stochastic_arrow import reenact_events

history = reenact_events(stoichiometric_matrix, result['events'], state)
```

## Building

    > make clean compile

This builds the extension package and installs it in editable mode.

**NOTE:** `make compile` without an explicit `clean` might not fully build the extension.

## Testing

`stochastic_arrow` uses [pytest](https://docs.pytest.org/en/latest/).
To run the main tests, in the source tree:

    > make test

or

    > pytest

There are additional command line features in test_arrow:

    > python -m test.test_arrow --help
    > python -m test.test_arrow --complexation
    > python -m test.test_arrow --complexation --runs 3
    > python -m test.test_arrow --obsidian
    > python -m test.test_arrow --memory
    > python -m test.test_arrow --time
    > python -m test.test_arrow --pickle
    > python -m test.test_arrow --test-fail-flagella
    > python -m test.test_arrow --test-fail-stdout
    > python -m test.test_hang

This test requires installing a version of matplotlib that's compatible with the installed numpy:

    > python -m test.test_arrow --plot

More examples:

    > pytest -k flagella

## Changelog

### Version 1.1.0
* Update build toolchain and automatically build/publish wheels for all
major platforms and recent Python versions.
* Build wheels with Numpy 2+ support

### Version 1.0.0

* Rename module to `stochastic_arrow` to avoid name conflict (Issue #51). **All users must update their import statements to use `stochastic_arrow` instead of `arrow`.**

### Version 0.5.2

* Update to Cython 0.29.34. (Cython 3.0.0 is now in beta.)

### Version 0.5.1

* Update to Cython 3.0.0a11 for compatibility with Python 3.11.
  Add `arrow.pxd` to work around a Cython 3.0.0 bug.
* Stop using deprecated `numpy.distutils` to avoid warnings and prepare for its
  removal in Python 3.12.
* Make `test_arrow.py --plot` compatible with Python 3.
* Fix `PytestReturnNotNoneWarning` warnings from pytest 7.2.0.

### Version 0.5.0

* Add the arrow_hang unit test which catches a nasty edge-case (Issue #48),
  fix the bug, and make the code more robust to some other potential bugs.

### Version 0.4.4

* Can pickle StochasticSystem instances.

### Version 0.3.0

* Introduced backwards-incompatible API change for supplying rates at `evolve()` time rather than `__init__()` for `StochasticSystem`.
