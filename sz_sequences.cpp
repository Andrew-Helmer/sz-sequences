/*
 * Copyright (C) Andrew Helmer 2026.
 * Licensed under MIT No Attribution (MIT-0) License: see LICENSE.
 *
 * If you use this code, I'd ask that you cite this GitHub repository in your
 * source code: https://github.com/Andrew-Helmer/sz-sequences
 *
 * This file implements the SZ sequence from "SZ Sequences: Binary-Based
 * (0, 2^q)-Sequences" by Ahmed et al. (2025), with shuffling and scrambling
 * from "Practical Hash-based Owen Scrambling" by Burley (2020).
 *
 * The generator matrices are computed at compile-time, using the direct
 * construction of the SZ sequence matrices from
 * https://andrew-helmer.github.io/sz-matrix-construction/
 * Nothing here allocates heap memory.
 */
#include "sz_sequences.h"

#include <cassert>
#include <cstdint>

namespace sz_sequences {

// The number of dimensions of the SZ sequence that can be generated. With
// 32x32 matrices, only the first 65536 dimensions are unique (after that, the
// generator matrices repeat). Higher values will take more time (and memory) at
// compile-time. With clang's default constexpr limits, up to 512 dimensions
// will work, but higher values may need a higher limit (e.g. -fconstexpr-steps
// with clang).
static constexpr int kNumDimensions = 16;
static_assert(kNumDimensions >= 2 && kNumDimensions <= 65536,
              "kNumDimensions must be between 2 and 65536.");

// The generator matrices are 32x32, so every coordinate has 32 bits.
static constexpr int kMatrixSize = 32;

/*
 * From here on down are the constexpr functions to construct the generator
 * matrices at compile-time, from "A Direct Construction of SZ Sequence
 * Matrices" (https://andrew-helmer.github.io/sz-matrix-construction/).
 */

// Gets the number of bits needed to represent a value, e.g. 15 -> 4 and
// 16 -> 5. For d > 0, this is floor(log_2(d)) + 1, where j = floor(log_2(d))
// is the index of the highest set bit of d (Step 3 of the algorithm).
static constexpr int GetBitLength(uint32_t value) {
  int bit_length = 0;
  while (value != 0) {
    bit_length++;
    value >>= 1;
  }
  return bit_length;
}

// To generate kNumDimensions (up to dimension d = kNumDimensions - 1), we need
// j + 1 basis matrices E_0, ..., E_j, where j = floor(log_2(d)). For 16
// dimensions (max d = 15), j = 3, so we need 4 basis matrices.
static constexpr int kNumBasisMatrices = GetBitLength(kNumDimensions - 1);

// A square column-major binary matrix of up to 32x32 entries. During matrix
// construction, each column is stored as a single 32-bit integer, where the
// entry in row r is bit r. Smaller matrices (like the basis matrices) just use
// the first columns and the lowest bits.
struct BinaryMatrix {
  uint32_t cols[kMatrixSize] = {};
};

static constexpr BinaryMatrix GetIdentityMatrix(const int size) {
  BinaryMatrix identity;
  for (int col = 0; col < size; col++) {
    identity.cols[col] = 1u << col;
  }
  return identity;
}

// Repeats a block matrix along the main diagonal, up to a target size. All
// other entries are zero. Used in Step 2.2.1 of the basis matrix algorithm
// to construct E'_i from E_i.
static constexpr BinaryMatrix RepeatDiagonal(const BinaryMatrix& block,
                                             const int block_size,
                                             const int target_size) {
  BinaryMatrix repeated;
  for (int offset = 0; offset < target_size; offset += block_size) {
    for (int col = 0; col < block_size; col++) {
      repeated.cols[offset + col] = block.cols[col] << offset;
    }
  }
  return repeated;
}

// Multiplies two size x size binary matrices, with all values modulo 2.
//
// In standard matrix multiplication, entry (k, col) of left * right is the sum
// over row of left[k][row] * right[row][col]. Viewed a whole column at a time,
// column "col" of the product is the sum of left.cols[row] for every row where
// entry (row, col) of the right matrix is 1. Since addition modulo 2 is xor,
// we just xor those columns together.
static constexpr BinaryMatrix MultiplyMatrices(const BinaryMatrix& left,
                                               const BinaryMatrix& right,
                                               const int size) {
  BinaryMatrix product;
  for (int col = 0; col < size; col++) {
    for (int row = 0; row < size; row++) {
      if ((right.cols[col] >> row) & 1) {
        product.cols[col] ^= left.cols[row];
      }
    }
  }
  return product;
}

// The basis matrices (E_0, E_1, ..., E_j) are the building blocks of all the
// generator matrices. The first one is the 1x1 matrix [1], and each
// power-of-two group of basis matrices doubles in size: 1x1, 2x2, 4x4, 4x4,
// 8x8, 8x8, 8x8, 8x8, etc.
//
// Mathematically, each doubling step builds a quadratic finite-field extension
// over GF(2) in matrix form (GF(2) -> GF(4) -> GF(16) -> GF(256) ...), using
// the iterated quadratic extension tower from "An Iterated Quadratic Extension
// of GF(2)" by Wiedemann (1988).
struct BasisMatrices {
  BinaryMatrix matrices[kNumBasisMatrices];
  int sizes[kNumBasisMatrices] = {};
};

// Generates the j + 1 basis block matrices E_0, E_1, ..., E_j, following the
// basis matrix algorithm in "New Matrix Construction Algorithm" in the blog
// post.
static constexpr BasisMatrices GenerateBasisMatrices() {
  BasisMatrices basis;
  // Step 1: Initialize E_0 = [1] and n = 1 (here n is subblock_size).
  basis.matrices[0] = GetIdentityMatrix(1);
  basis.sizes[0] = 1;

  // Step 2: While n <= j (with Step 2.3, n <- 2n, at the end of each loop).
  // In each iteration, subblock_size (n) is both the number of basis matrices
  // generated so far, and the size of the last basis matrix (beta).
  for (int subblock_size = 1; subblock_size < kNumBasisMatrices;
       subblock_size *= 2) {
    // Step 2.1: Set beta = E_{n-1} (which is an n x n matrix).
    const BinaryMatrix beta = basis.matrices[subblock_size - 1];

    // Step 2.2: For i = 0 ... n - 1 (while n + i <= j):
    for (int i = 0;
         i < subblock_size && subblock_size + i < kNumBasisMatrices;
         i++) {
      // Step 2.2.1: Let E'_i ("enlarged") be E_i repeated along the main
      // diagonal to size n x n.
      const BinaryMatrix enlarged = RepeatDiagonal(
          basis.matrices[i], basis.sizes[i], subblock_size);
      const BinaryMatrix beta_enlarged =
          MultiplyMatrices(beta, enlarged, subblock_size);

      // Step 2.2.2: Construct the 2n x 2n basis matrix E_{n+i} as the 2 x 2
      // block matrix:
      //   E_{n+i} = | 0_n        E'_i |
      //             | beta*E'_i  E'_i |
      // (This comes from multiplying the companion matrix of x^2 + x + beta = 0
      // by E'_i repeated on the diagonal.)
      BinaryMatrix& next_basis = basis.matrices[subblock_size + i];
      for (int col = 0; col < subblock_size; col++) {
        next_basis.cols[col] = beta_enlarged.cols[col] << subblock_size;
        next_basis.cols[subblock_size + col] =
            enlarged.cols[col] | (enlarged.cols[col] << subblock_size);
      }
      basis.sizes[subblock_size + i] = subblock_size * 2;
    }
  }
  return basis;
}

// Constructs the 32x32 upper triangular block-Pascal matrix P(X) from an
// m x m square block X ("block" of size "block_size"), where each m x m block
// in block row i and block column j (for j >= i) is defined as:
//   P(X)_{i,j} = binom(j, i) * X^(j - i) (mod 2).
static constexpr BinaryMatrix GetBlockPascalMatrix(const BinaryMatrix& block,
                                                   const int block_size) {
  const int num_blocks = kMatrixSize / block_size;
  BinaryMatrix pascal_matrix;

  // We loop over each block diagonal (power = j - i), which all have the same
  // block power X^(j - i).
  BinaryMatrix block_power = GetIdentityMatrix(block_size);
  for (int power = 0; power < num_blocks; power++) {
    for (int block_row = 0; block_row + power < num_blocks; block_row++) {
      const int block_col = block_row + power;
      // We need to know if binom(j, i) = binom(block_col, block_row) is odd
      // (1 mod 2). Lucas' theorem, from "Theorie des Fonctions Numeriques
      // Simplement Periodiques" by Lucas (1878), says that binom(n, k) mod 2
      // is the product of binom(n_b, k_b) over each bit b of n and k. For
      // single bits, binom(0, 1) = 0 and all other cases are 1, so binom(n, k)
      // is odd if and only if every set bit in k is also set in n.
      if ((block_col & block_row) != block_row) {
        continue;
      }

      for (int col = 0; col < block_size; col++) {
        pascal_matrix.cols[block_col * block_size + col] |=
            block_power.cols[col] << (block_row * block_size);
      }
    }
    block_power = MultiplyMatrices(block_power, block, block_size);
  }
  return pascal_matrix;
}

// Multiplies two size x size non-singular upper-triangular (NUT) binary
// matrices modulo 2. This is the same as MultiplyMatrices, except that because
// the right matrix is upper-triangular, entry (row, col) can only be 1 when
// row <= col. Skipping the lower triangle cuts the work in half.
static constexpr BinaryMatrix MultiplyNUTMatrices(const BinaryMatrix& left,
                                                  const BinaryMatrix& right,
                                                  const int size) {
  BinaryMatrix product;
  for (int col = 0; col < size; col++) {
    for (int row = 0; row <= col; row++) {
      if ((right.cols[col] >> row) & 1) {
        product.cols[col] ^= left.cols[row];
      }
    }
  }
  return product;
}

static constexpr uint32_t ReverseBits(uint32_t value) {
  value = (((value & 0xaaaaaaaa) >> 1) | ((value & 0x55555555) << 1));
  value = (((value & 0xcccccccc) >> 2) | ((value & 0x33333333) << 2));
  value = (((value & 0xf0f0f0f0) >> 4) | ((value & 0x0f0f0f0f) << 4));
  value = (((value & 0xff00ff00) >> 8) | ((value & 0x00ff00ff) << 8));
  return ((value >> 16) | (value << 16));
}

struct GeneratorMatrices {
  BinaryMatrix matrices[kNumDimensions];
};

// Generates all the SZ sequence generator matrices, using the algorithm from
// https://andrew-helmer.github.io/sz-matrix-construction/
//
// In the single-matrix algorithm from the blog post, the N x N generator
// matrix M for dimension d is constructed by:
//   Step 1: Initialize M = I_N (and if d = 0, return M in Step 2).
//   Step 3-4: Generate the basis block matrices E_0, E_1, ..., E_j, where
//             j = floor(log_2(d)).
//   Step 5: For i = 0 ... j, if bit i of d is set ((d >> i) & 1), multiply M
//           on the right by the block-Pascal matrix P(E_i): M <- M * P(E_i).
//
// As described in the "Generating Multiple Matrices" section of the blog post,
// we can generate all dimensions at once without repeating work: for each i
// from 0 to j, dimensions d = 2^i + prev_d have the same set bits as prev_d
// plus bit i, so their generator matrices are just M_{prev_d} * P(E_i). Each
// P(E_i) therefore doubles the array of generator matrices by right-multiplying
// the previously generated ones.
static constexpr GeneratorMatrices GenerateSZMatrices() {
  GeneratorMatrices sz;
  // Step 1-2: Initialize M_0 = I_N (the 32x32 identity matrix).
  sz.matrices[0] = GetIdentityMatrix(kMatrixSize);

  // Step 4: Generate the j + 1 basis block matrices E_0, ..., E_j.
  const BasisMatrices basis = GenerateBasisMatrices();

  // Step 5: For i = 0 ... j, compute P(E_i) and right-multiply each previously
  // generated matrix M by P(E_i) (M <- M * P(E_i)). Because both M and P(E_i)
  // are non-singular upper-triangular (NUT) matrices, we can use
  // MultiplyNUTMatrices.
  int num_generated = 1;
  for (int i = 0; i < kNumBasisMatrices; i++) {
    const BinaryMatrix pascal_matrix =
        GetBlockPascalMatrix(basis.matrices[i], basis.sizes[i]);

    const int prev_len = num_generated;
    for (int prev_d = 0;
         prev_d < prev_len && num_generated < kNumDimensions; prev_d++) {
      sz.matrices[num_generated] = MultiplyNUTMatrices(
          sz.matrices[prev_d], pascal_matrix, kMatrixSize);
      num_generated++;
    }
  }

  // Reverse the bits of each column so that the first row (row 0) gives the
  // most significant bit (bit 31) of the 32-bit coordinate.
  for (int dim = 0; dim < kNumDimensions; dim++) {
    for (int col = 0; col < kMatrixSize; col++) {
      sz.matrices[dim].cols[col] = ReverseBits(sz.matrices[dim].cols[col]);
    }
  }
  return sz;
}

static constexpr GeneratorMatrices kSZMatrices = GenerateSZMatrices();

// Gets a coordinate of the SZ sequence as a 32-bit fixed point value.
static inline uint32_t GetSZBits(const uint32_t index, const int dim) {
  assert(dim >= 0 && dim < kNumDimensions);
  // In a digital sequence, the coordinate is generated by multiplying the
  // generator matrix by the column vector of bits of the sample index. That
  // matrix-vector product is the xor of the generator matrix columns for each
  // set bit of the index.
  uint32_t bits = 0;
  uint32_t remaining_index = index;
  for (int bit = 0; remaining_index != 0; bit++) {
    if (remaining_index & 1) {
      bits ^= kSZMatrices.matrices[dim].cols[bit];
    }
    remaining_index >>= 1;
  }
  return bits;
}

// Converts a 32-bit fixed point value to a float in [0, 1). A float only has 24
// bits of precision, so we only use the 24 most significant bits. Otherwise
// values very close to 1 could be rounded up to exactly 1.
static inline float BitsToFloat(const uint32_t bits) {
  return static_cast<float>(bits >> 8) / 16777216.0f;  // 2^24
}

float GetSZCoord(const int index, const int dim) {
  const uint32_t unscrambled_bits =
      GetSZBits(static_cast<uint32_t>(index), dim);
  return BitsToFloat(unscrambled_bits);
}

/*
 * Shuffling and scrambling functions from "Practical Hash-based Owen
 * Scrambling" (Burley 2020).
 */

static inline uint32_t Hash(uint32_t value) {
  // finalizer from murmurhash3
  value ^= value >> 16;
  value *= 0x85ebca6bu;
  value ^= value >> 13;
  value *= 0xc2b2ae35u;
  value ^= value >> 16;
  return value;
}

// Burley's version of the permutation from "Stratified Sampling for Stochastic
// Transparency" (Laine and Karras 2011). Because each even constant is
// multiplied by value, carries and bit-mixes only flow from lower bits into
// higher bits, so each bit of the output only depends on the seed and the same
// bit or lower bits of the input.
static inline uint32_t LaineKarrasPermutation(uint32_t value,
                                              const uint32_t seed) {
  value += seed;
  value ^= value * 0x6c50b47cu;
  value ^= value * 0xb82f1e52u;
  value ^= value * 0xc7afe638u;
  value ^= value * 0x8d22f6e6u;
  return value;
}

// Owen-scrambling in base-2, from "Randomly Permuted (t,m,s)-Nets and
// (t,s)-Sequences" (Owen 1995). In base-2 Owen-scrambling, whether a bit is
// flipped depends on all the more significant bits above it. By reversing the
// bits before and after LaineKarrasPermutation, each bit of the output only
// depends on the seed and the same bit or higher bits of the input.
static inline uint32_t NestedUniformScrambleBase2(uint32_t value,
                                                  const uint32_t seed) {
  value = ReverseBits(value);
  value = LaineKarrasPermutation(value, seed);
  value = ReverseBits(value);
  return value;
}

static inline uint32_t CombineHashes(const uint32_t seed,
                                     const uint32_t value) {
  return seed ^ (value + (seed << 6) + (seed >> 2));
}

float GetShuffledScrambledSZCoord(const int index,
                                  const int dim,
                                  const uint32_t seed) {
  // Hashing the seed first means that similar seeds (e.g. 0, 1, 2...) will
  // still give very different shuffling and scrambling.
  const uint32_t hashed_seed = Hash(seed);

  // Shuffling is Owen-scrambling on the index. Any power-of-two disjoint
  // subsequence maps to another power-of-two disjoint subsequence, so the
  // stratification is maintained. The shuffling must be the same for all
  // dimensions, so it doesn't depend on dim.
  const uint32_t shuffled_index =
      NestedUniformScrambleBase2(static_cast<uint32_t>(index), hashed_seed);

  // Scrambling is Owen-scrambling on the coordinate, with a different seed for
  // each dimension.
  const uint32_t dim_seed =
      CombineHashes(hashed_seed, Hash(static_cast<uint32_t>(dim)));
  const uint32_t unscrambled_bits = GetSZBits(shuffled_index, dim);
  const uint32_t scrambled_bits =
      NestedUniformScrambleBase2(unscrambled_bits, dim_seed);
  return BitsToFloat(scrambled_bits);
}

}  // namespace sz_sequences
