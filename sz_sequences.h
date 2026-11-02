/*
 * Copyright (C) Andrew Helmer 2026.
 * Licensed under MIT No Attribution (MIT-0) License: see LICENSE.
 *
 * If you use this code, I'd ask that you cite this GitHub repository in your
 * source code: https://github.com/Andrew-Helmer/sz-sequences
 *
 * This header defines functions to get points from the SZ sequence, from
 * "SZ Sequences: Binary-Based (0, 2^q)-Sequences" by Ahmed et al. (2025),
 * either unscrambled and unshuffled, or shuffled and scrambled using the
 * hash-based approach from "Practical Hash-based Owen Scrambling"
 * by Burley (2020).
 *
 * The SZ sequence is a base-2 sequence where every aligned pair of dimensions
 * (0+1, 2+3, etc.) forms a base-2 (0,2)-sequence, every aligned group of four
 * dimensions (0-3, 4-7, etc.) forms a base-4 (0,4)-sequence, every aligned
 * group of 16 dimensions forms a base-16 (0,16)-sequence, and so on.
 *
 * The number of dimensions that can be generated is set by kNumDimensions in
 * sz_sequences.cpp, which is 16 by default. All dim values must be less than
 * that.
 */
#ifndef SZ_SEQUENCES_H
#define SZ_SEQUENCES_H

#include <cstdint>

namespace sz_sequences {

// Gets a single coordinate of the unscrambled, unshuffled SZ sequence, in
// [0, 1).
float GetSZCoord(const int index, const int dim);

// Gets a single coordinate of the shuffled and scrambled SZ sequence, in
// [0, 1). The same seed must be used for all dimensions of a point, otherwise
// the dimensions will be shuffled differently and the points won't be
// stratified across dimensions. Different seeds will decorrelate the
// sequences.
float GetShuffledScrambledSZCoord(const int index,
                                  const int dim,
                                  const uint32_t seed);

}  // namespace sz_sequences

#endif  // SZ_SEQUENCES_H
