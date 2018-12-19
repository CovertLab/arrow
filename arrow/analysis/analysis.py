
from __future__ import absolute_import, division, print_function

import numpy as np

def moving_average(time, counts, sampling_time, scale):
	# TODO: vectorize

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
