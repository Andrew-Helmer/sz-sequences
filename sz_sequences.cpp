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
 * The generator matrices are computed at compile-time, using a new direct
 * construction of the SZ sequence matrices. Nothing here allocates heap memory.
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

namespace {

/*
 * From here on down are the constexpr functions to construct the generator
 * matrices at compile-time.
 */

// Gets the number of bits needed to represent a value, e.g. 15 -> 4 and
// 16 -> 5.
constexpr int GetBitLength(uint32_t value) {
  int bit_length = 0;
  while (value != 0) {
    bit_length++;
    value >>= 1;
  }
  return bit_length;
}

// Each power-of-two number of dimensions needs one more basis matrix, so for
// 16 dimensions we need 4 basis matrices.
static constexpr int kNumBasisMatrices = GetBitLength(kNumDimensions - 1);

// A square binary matrix of up to 32x32 entries. During matrix construction,
// each row is stored as a single 32-bit integer, where the entry in column c
// is bit c. Smaller matrices (like the basis matrices) just use the first rows
// and the lowest bits.
struct BinaryMatrix {
  uint32_t rows[kMatrixSize] = {};
};

constexpr BinaryMatrix GetIdentityMatrix(const int size) {
  BinaryMatrix identity;
  for (int row = 0; row < size; row++) {
    identity.rows[row] = 1u << row;
  }
  return identity;
}

// Multiplies two size x size binary matrices, with all values modulo 2.
//
// In standard matrix multiplication, entry (row, k) of left * right is the sum
// over col of left[row][col] * right[col][k]. Viewed a whole row at a time,
// row "row" of the product is the sum of right.rows[col] for every col where
// entry (row, col) of the left matrix is 1. Since addition modulo 2 is xor, we
// just xor those rows together.
constexpr BinaryMatrix MultiplyMatrices(const BinaryMatrix& left,
                                        const BinaryMatrix& right,
                                        const int size) {
  BinaryMatrix product;
  for (int row = 0; row < size; row++) {
    for (int col = 0; col < size; col++) {
      if ((left.rows[row] >> col) & 1) {
        product.rows[row] ^= right.rows[col];
      }
    }
  }
  return product;
}

// Repeats a block matrix along the main diagonal, up to a target size. All
// other entries are zero.
constexpr BinaryMatrix RepeatDiagonal(const BinaryMatrix& block,
                                      const int block_size,
                                      const int target_size) {
  BinaryMatrix repeated;
  for (int offset = 0; offset < target_size; offset += block_size) {
    for (int row = 0; row < block_size; row++) {
      repeated.rows[offset + row] = block.rows[row] << offset;
    }
  }
  return repeated;
}

// The basis matrices are the building blocks of all the generator matrices.
// The first one is the 1x1 matrix [1], and each power-of-two group of basis
// matrices doubles in size: 1x1, 2x2, 4x4, 4x4, 8x8, 8x8, 8x8, 8x8, etc.
//
// Mathematically, each doubling step builds a quadratic finite-field extension
// over GF(2) in matrix form (GF(2) -> GF(4) -> GF(16) -> GF(256) ...), using
// the iterated quadratic extension tower from "An Iterated Quadratic Extension
// of GF(2)" by Wiedemann (1988).
struct BasisMatrices {
  BinaryMatrix matrices[kNumBasisMatrices];
  int sizes[kNumBasisMatrices] = {};
};

constexpr BasisMatrices GenerateBasisMatrices() {
  BasisMatrices basis;
  basis.matrices[0] = GetIdentityMatrix(1);
  basis.sizes[0] = 1;

  // In each iteration, subblock_size is both the number of basis matrices
  // generated so far, and the size of the last basis matrix (beta).
  for (int subblock_size = 1; subblock_size < kNumBasisMatrices;
       subblock_size *= 2) {
    const BinaryMatrix beta = basis.matrices[subblock_size - 1];
    for (int idx = 0;
         idx < subblock_size && subblock_size + idx < kNumBasisMatrices;
         idx++) {
      // Each of the earlier basis matrices is repeated along the diagonal to
      // subblock_size, then multiplied by the companion matrix of
      // x^2 + x + beta = 0 to construct a new basis matrix of twice the size:
      //   | 0    I | * | enlarged 0        | = | 0             enlarged |
      //   | beta I |   | 0        enlarged |   | beta*enlarged enlarged |
      const BinaryMatrix enlarged = RepeatDiagonal(
          basis.matrices[idx], basis.sizes[idx], subblock_size);
      const BinaryMatrix beta_enlarged =
          MultiplyMatrices(beta, enlarged, subblock_size);

      BinaryMatrix& next_basis = basis.matrices[subblock_size + idx];
      for (int row = 0; row < subblock_size; row++) {
        next_basis.rows[row] = enlarged.rows[row] << subblock_size;
        next_basis.rows[subblock_size + row] =
            beta_enlarged.rows[row] | (enlarged.rows[row] << subblock_size);
      }
      basis.sizes[subblock_size + idx] = subblock_size * 2;
    }
  }
  return basis;
}

// Constructs the 32x32 upper triangular block-Pascal matrix from a square
// block, where the block in block row i and block column j (for j >= i) is
// binom(j, i) * block^(j - i), modulo 2.
constexpr BinaryMatrix GetBlockPascalMatrix(const BinaryMatrix& block,
                                            const int block_size) {
  const int num_blocks = kMatrixSize / block_size;
  BinaryMatrix pascal_matrix;

  // We loop over each block diagonal, which all have the same block power.
  BinaryMatrix block_power = GetIdentityMatrix(block_size);
  for (int power = 0; power < num_blocks; power++) {
    for (int block_row = 0; block_row + power < num_blocks; block_row++) {
      const int block_col = block_row + power;
      // We need to know if binom(block_col, block_row) is odd (1 mod 2).
      // Lucas' theorem, from "Theorie des Fonctions Numeriques Simplement
      // Periodiques" by Lucas (1878), says that binom(n, k) mod 2 is the
      // product of binom(n_i, k_i) over each bit i of n and k. For single
      // bits, binom(0, 1) = 0 and all other cases are 1, so binom(n, k) is
      // odd if and only if every set bit in k is also set in n.
      if ((block_col & block_row) != block_row) {
        continue;
      }

      for (int row = 0; row < block_size; row++) {
        pascal_matrix.rows[block_row * block_size + row] |=
            block_power.rows[row] << (block_col * block_size);
      }
    }
    block_power = MultiplyMatrices(block_power, block, block_size);
  }
  return pascal_matrix;
}

// For fast point generation, each 32x32 generator matrix is stored as an array
// of 32 column integers, where row r of a column is stored in bit (31 - r), so
// that the first row (row 0) gives the most significant bit of the coordinate.
struct GeneratorMatrices {
  uint32_t columns[kNumDimensions][kMatrixSize] = {};
};

// Generates all the SZ sequence generator matrices. The first matrix is the
// identity, and each block-Pascal matrix (one for each basis matrix) doubles
// the number of generator matrices by right-multiplying all of the previous
// generator matrices.
//
// Because block-Pascal matrices satisfy P(A + B) = P(A) * P(B) for commuting
// matrices A and B modulo 2, right-multiplying by each basis block-Pascal
// matrix generates the block-Pascal matrices for all combinations of the basis
// matrices.
constexpr GeneratorMatrices GenerateSZMatrices() {
  GeneratorMatrices sz;
  for (int col = 0; col < kMatrixSize; col++) {
    sz.columns[0][col] = 0x80000000u >> col;
  }

  const BasisMatrices basis = GenerateBasisMatrices();
  int num_generated = 1;
  for (int basis_idx = 0; basis_idx < kNumBasisMatrices; basis_idx++) {
    const BinaryMatrix pascal_matrix =
        GetBlockPascalMatrix(basis.matrices[basis_idx], basis.sizes[basis_idx]);

    const int prev_len = num_generated;
    for (int prev_idx = 0;
         prev_idx < prev_len && num_generated < kNumDimensions; prev_idx++) {
      // Just as multiplying on the left combines the rows of the right matrix,
      // multiplying on the right (prev_matrix * pascal_matrix) combines the
      // columns of the left matrix: column "col" of the product is the sum
      // (xor) of prev_cols[row] for every row where entry (row, col) of the
      // Pascal matrix is 1. Because the Pascal matrix is upper triangular,
      // entry (row, col) can only be 1 when row <= col.
      const uint32_t* prev_cols = sz.columns[prev_idx];
      uint32_t* next_cols = sz.columns[num_generated];
      for (int col = 0; col < kMatrixSize; col++) {
        for (int row = 0; row <= col; row++) {
          if ((pascal_matrix.rows[row] >> col) & 1) {
            next_cols[col] ^= prev_cols[row];
          }
        }
      }
      num_generated++;
    }
  }
  return sz;
}

static constexpr GeneratorMatrices kSZMatrices = GenerateSZMatrices();

/*
 * Shuffling and scrambling functions from "Practical Hash-based Owen
 * Scrambling" (Burley 2020).
 */
inline uint32_t ReverseBits(uint32_t value) {
  value = (((value & 0xaaaaaaaa) >> 1) | ((value & 0x55555555) << 1));
  value = (((value & 0xcccccccc) >> 2) | ((value & 0x33333333) << 2));
  value = (((value & 0xf0f0f0f0) >> 4) | ((value & 0x0f0f0f0f) << 4));
  value = (((value & 0xff00ff00) >> 8) | ((value & 0x00ff00ff) << 8));
  return ((value >> 16) | (value << 16));
}

// Burley's version of the permutation from "Stratified Sampling for Stochastic
// Transparency" (Laine and Karras 2011). Because each even constant is
// multiplied by value, carries and bit-mixes only flow from lower bits into
// higher bits, so each bit of the output only depends on the seed and the same
// bit or lower bits of the input.
inline uint32_t LaineKarrasPermutation(uint32_t value, const uint32_t seed) {
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
inline uint32_t NestedUniformScrambleBase2(uint32_t value,
                                           const uint32_t seed) {
  value = ReverseBits(value);
  value = LaineKarrasPermutation(value, seed);
  value = ReverseBits(value);
  return value;
}

inline uint32_t CombineHashes(const uint32_t seed, const uint32_t value) {
  return seed ^ (value + (seed << 6) + (seed >> 2));
}

inline uint32_t Hash(uint32_t value) {
  // finalizer from murmurhash3
  value ^= value >> 16;
  value *= 0x85ebca6bu;
  value ^= value >> 13;
  value *= 0xc2b2ae35u;
  value ^= value >> 16;
  return value;
}

// Converts a 32-bit fixed point value to a float in [0, 1). A float only has 24
// bits of precision, so we only use the 24 most significant bits. Otherwise
// values very close to 1 could be rounded up to exactly 1.
inline float BitsToFloat(const uint32_t bits) {
  return static_cast<float>(bits >> 8) / 16777216.0f;  // 2^24
}

// Gets a coordinate of the SZ sequence as a 32-bit fixed point value.
inline uint32_t GetSZBits(const uint32_t index, const int dim) {
  assert(dim >= 0 && dim < kNumDimensions);
  // In a digital sequence, the coordinate is generated by multiplying the
  // generator matrix by the column vector of bits of the sample index. That
  // matrix-vector product is the xor of the generator matrix columns for each
  // set bit of the index.
  uint32_t bits = 0;
  uint32_t remaining_index = index;
  for (int bit = 0; remaining_index != 0; bit++) {
    if (remaining_index & 1) {
      bits ^= kSZMatrices.columns[dim][bit];
    }
    remaining_index >>= 1;
  }
  return bits;
}
}  // namespace

float GetSZCoord(const int index, const int dim) {
  return BitsToFloat(GetSZBits(static_cast<uint32_t>(index), dim));
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
  const uint32_t scrambled_bits =
      NestedUniformScrambleBase2(GetSZBits(shuffled_index, dim), dim_seed);
  return BitsToFloat(scrambled_bits);
}

}  // namespace sz_sequences
