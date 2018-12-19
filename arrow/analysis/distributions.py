'''

This submodule is for statistical reductions on stochastic simulation output.

'''

from __future__ import absolute_import, division, print_function

import numpy as np

def moving_average(time, counts, sampling_time, scale):
	'''
	Compute a moving average on stochastic simulation output.

	Parameters
	----------

	time: 1D ndarray of floats
		The times at which each timepoint occurs.  Assumed to be strictly
		increasing.

	counts: 2D ndarray of integers
		Timepoints should be along the first axis (rows) while individual
		molecular species should be along the second axis (columns).

	sampling_time: 1D ndarray of floats
		The new times at which to sample while generating the moving average.
		Again, assumed to be strictly increasing.

	scale: positive float
		The time scale of the moving average.  Larger values result in
		'smoother' curves.

	Returns
	-------

	averaged: 2D ndarray of floats
		The moving average of the counts.

	Notes
	-----

	This function was implemented specifically for stochastic simulation output
	but is likely applicable to any data.

	As implemented, this moving average will always lag behind the 'true'
	values.  I.e. an oscillatory trajectory will exhibit a phase shift.

	Too small a scale will appear no different from the unaveraged counts,
	while too long a scale will lead to extreme lag and wash out any critical
	magnitude variations.

	This could be sped up considerably if compiled via JIT, or expressed using
	other compilable operations (e.g. a convolution).
	'''

	counts = np.asarray(counts)

	n_sampling = sampling_time.size

	if counts.ndim == 1:
		shape = (n_sampling,)

	else:
		shape = (n_sampling, counts.shape[1])

	averaged = np.empty(shape, np.float64)

	for i in xrange(n_sampling):
		j = _last_where(sampling_time[i] >= time)

		if i == 0:
			v = counts[j]

		else:
			dt = sampling_time[i] - sampling_time[i-1]

			w = dt / scale

			v = (w*counts[j] + averaged[i-1]) / (w+1)

		averaged[i] = v

	return averaged

def _last_where(bool_array):
	'''
	Return the index of the last element that is true in a 1D boolean array.

	This could be made much faster and more robust.
	'''
	return np.where(bool_array)[0][-1]

if __name__ == '__main__':
	import matplotlib.pyplot as plt

	N = 300
	Np = 3000

	t = np.random.exponential(size = N).cumsum()
	c = (2*(np.random.random(size = N) > 0.5)-1).cumsum()+1000

	tp = np.linspace(t[0], t[-1], Np)

	for scale in (1.0, 5.0, 25.0, 125.0):
		cp = moving_average(t, c, tp, scale)
		plt.plot(tp, cp, lw = 2.0, label = '{:0.2f}'.format(scale))

	plt.step(
		t, c,
		where = 'post',
		c = 'k', lw = 1.0, #alpha = 0.5,
		label = 'exact'
		)

	plt.legend(loc = 'best')

	plt.show()
