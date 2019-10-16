from __future__ import absolute_import, division, print_function

from six import moves

def choose(n, k):
    '''
    Enumerate the number of ways to pick 'k' items from a set of 'n', assuming
    that items are unique but without concern for order.  E.g. (A, B) and
    (B, A) from {A, B, C} are only counted once.
    '''

    combinations = 1.0

    for i in moves.range(k):
        combinations *= (n-i)/(i+1)

    return combinations

def multichoose(n, k):
    '''
    Enumerate the number of ways to pick elements from different sets of items,
    where 'n' is a vector of the number of items in each set, and 'k' is a
    vector of the number of items to choose from each set.
    '''

    assert n.size == k.size

    combinations = 1.0

    for i in moves.range(n.size):
        combinations *= choose(n[i], k[i])

    return combinations
