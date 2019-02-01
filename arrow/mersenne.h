/*
 * The Mersenne Twister pseudo-random number generator (PRNG)
 *
 * This is an implementation of fast PRNG called MT19937, meaning it has a
 * period of 2^19937-1, which is a Mersenne prime.
 *
 * This PRNG is fast and suitable for non-cryptographic code.  For instance, it
 * would be perfect for Monte Carlo simulations, etc.
 *
 * For all the details on this algorithm, see the original paper:
 * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/ARTICLES/mt.pdf
 *
 * Written by Christian Stigen Larsen
 * Distributed under the modified BSD license.
 * 2015-02-17, 2017-12-06
 */

#ifndef MERSENNE_H
#define MERSENNE_H

#include <stdint.h>

// static const size_t TWISTER_SIZE = 624;
#define TWISTER_SIZE 624

typedef struct MTState MTState;
struct MTState {
  uint32_t MT[TWISTER_SIZE];
  uint32_t MT_TEMPERED[TWISTER_SIZE];
  size_t index;
};

// Initialize Mersenne Twister with given seed value.
void seed(MTState *state, uint32_t seed_value);

// Extract a pseudo-random unsigned 32-bit integer in the range 0 ... UINT32_MAX
uint32_t rand_u32(MTState *state);

// Sample from a uniform distribution
double sample_uniform(MTState *state);

// Sample from an exponential distribution of a given lambda
double sample_exponential(MTState *state, double lambda);

#endif
