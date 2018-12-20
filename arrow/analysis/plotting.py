'''

This submodule is for plotting functions; i.e., anything dependent upon
matplotlib.

'''

from __future__ import absolute_import, division, print_function

import matplotlib.pyplot as plt

def plot_full_history(axes, time, counts, *args, **kwargs):
    '''
    Light wrapper around matplotlib's 'step' plotting function.  Specifically
    (and correctly) requests post-point stepping, rather than pre-point
    stepping, which is the default, and would incorrectly suggest an
    instantaneous change in abundance at the beginning of the time step rather
    than at the end.
    '''
    axes.step(
        time,
        counts,
        *args,
        where = 'post',
        **kwargs
        )
