# SZ Sequences

This is a short and simple C++ implementation of the SZ sequence, from ["SZ Sequences: Binary-Based (0, 2^q)-Sequences"](https://dl.acm.org/doi/10.1145/3763272) by Ahmed et al. (2025). The points can be shuffled and scrambled using the hash-based approach from ["Practical Hash-based Owen Scrambling"](https://jcgt.org/published/0009/04/01/) by Burley (2020).

It's just one header and one cpp file. The generator matrices are all computed at compile-time using a new [direct construction](https://andrew-helmer.github.io/sz-matrix-construction/), and nothing allocates heap memory.

## The SZ Sequence

The SZ sequence is a base-2 sequence, like the Sobol' sequence, but with much better stratification across higher dimensions. Every aligned pair of dimensions (0+1, 2+3, etc.) forms a base-2 (0,2)-sequence, which means that they are stratified on all base-2 elementary intervals for a power-of-two prefix of samples. Additionally, every aligned group of four dimensions (0-3, 4-7, etc.) forms a base-4 (0,4)-sequence, and all 16 dimensions together form a base-16 (0,16)-sequence (and the sequence is infinitely extensible to higher dimensionality, e.g. 256 dimensions form a base-256 (0,256)-sequence).

So for example, the first 256 points will be stratified on the 16x16 grid in every pair of the 16 dimensions, and also on the 4x4x4x4 grid in each aligned group of four dimensions. Actually this is true for all power-of-two "disjoint subsequences". So it's not just the first 256 samples, but the next 256, and the 256 after that, and so on.

This makes the first four dimensions especially nice for rendering, where we often combine two inherently 2D sampling domains, like a 2D pixel position and a 2D lens position, or a 2D light sample and a 2D BRDF sample. Each pair gets the ideal 2D stratification, and the four dimensions together are stratified as well as possible.

## Using the Code

Just add sz_sequences.h and sz_sequences.cpp to your project. It needs C++14 or higher.

To get a 4D shuffled and scrambled SZ sequence point:

<pre><code>#include "sz_sequences.h"

float point[4];
for (int dim = 0; dim < 4; dim++) {
  point[dim] = sz_sequences::GetShuffledScrambledSZCoord(index, dim, seed);
}</code></pre>

Each coordinate will be in [0, 1). There is also a GetSZCoord function, to get the unscrambled, unshuffled points. This is mostly for academic purposes, I don't think there is a reason for a renderer to use this.

If you only want to use the 4D sequence with "padding" (different shufflings for different groups of 4D sequences), then your dimension should go 0,1,2,3,0,1,... repeating only the first four dimensions, and you should change the seed after every 4 dimensions.

### Performance Notes

Generating each coordinate individually involves some redundant work (shuffling the index, for example). So if it suits your renderer, you might want to combine all four samples into a single call.

It's also possible to significantly optimize Burley scrambling by precalculating the unscrambled bits for, say, 4K, 16K, or 64K, only shuffle the index within that range, and scramble the points. That avoids the bit reversal and matrix multiplications which are most (about 60%) of the CPU cost.

### Shuffling and Scrambling

Scrambling randomizes the positions of the points, and shuffling randomizes the order of the points, without breaking their progressive stratification properties. Both of these are Owen-scrambling in base 2, using the hash functions from Burley's paper. Shuffling is Owen-scrambling on the sample indices, and scrambling is Owen-scrambling on the coordinates.

Different seeds will decorrelate the sequences, so in a renderer you might use a different seed for each pixel. It's important to use the same seed for every dimension of a point, otherwise the dimensions will be shuffled differently and the points won't be stratified across dimensions.

Because of the shuffling, the first N points with a given seed won't be the first N points of the unshuffled sequence. But any power-of-two disjoint subsequence maps to another power-of-two disjoint subsequence, so it will still have all the same stratification properties.

### Number of Dimensions

The number of dimensions is set by kNumDimensions at the top of sz_sequences.cpp, which is 16 by default. The generator matrices take up 128 bytes per dimension (2KB for 16 dimensions). With 32x32 matrices, only the first 65536 dimensions are unique (after that, the generator matrices repeat).

Higher values will take more time to compute at compile-time. With clang's default constexpr limits, up to 512 dimensions will work, but higher values may need a higher limit (e.g. -fconstexpr-steps with clang).

## Licensing

Licensed under the MIT No Attribution (MIT-0) license. See the [LICENSE](/LICENSE) file. You can use this code for any purpose, and you don't need to include the license, but if you use it, I'd ask that you cite this GitHub repository in your source code: https://github.com/Andrew-Helmer/sz-sequences
