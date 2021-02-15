#ifndef SOBOLH
#define SOBOLH

#include "sobolmatrices.h"
#include "siphash.h"
#include <stdexcept>
#include <cmath>
static const float FloatOneMinusEpsilon = 0.99999994;

static inline uint32_t hash_combine(uint32_t seed, uint32_t v) {
  return seed ^ (v + (seed << 6) + (seed >> 2));
}

static inline  uint32_t reverse_bits(uint32_t x) {
  x = (((x & 0xaaaaaaaa) >> 1) | ((x & 0x55555555) << 1));
  x = (((x & 0xcccccccc) >> 2) | ((x & 0x33333333) << 2));
  x = (((x & 0xf0f0f0f0) >> 4) | ((x & 0x0f0f0f0f) << 4));
  x = (((x & 0xff00ff00) >> 8) | ((x & 0x00ff00ff) << 8));
  return ((x >> 16) | (x << 16));
}

namespace spacefillr {

static inline uint32_t hash_u32(uint32_t x, uint64_t seed1, uint64_t seed2) {
  uint64_t out;
  uint64_t k[] = {seed1, seed2};
  spacefillr::siphash((uint8_t *)(&x), 4, (uint8_t *)(k), (uint8_t *)(&out), 8);
  return out;
}

static inline uint32_t hash_u32(uint32_t n, uint32_t seed) {
  // Seeding.
  n = 0x6217c6e1 ^ (n + (seed * (0x9e3779b9)));

  // From https://github.com/skeeto/hash-prospector
  n ^= n >> 17;
  n = n * (0xed5ad4bb);
  n ^= n >> 11;
  n = n * (0xac4c1b51);
  n ^= n >> 15;
  n = n * (0x31848bab);
  n ^= n >> 14;
  return(n);
}

}//namespace spacefillr

//------------------------------------------------------

static inline  uint32_t owen_scramble_u32(uint32_t x, uint32_t seed) {
  uint32_t in_bits = x;
  uint32_t out_bits = x;

  // Do the Owen scramble.
  for (uint32_t bit = 0; bit < 31; ++bit) {
    uint32_t high_mask = ~((1 << (bit + 1)) - 1);
    uint32_t hash = spacefillr::hash_u32(in_bits & high_mask, seed, bit);
    out_bits ^= hash & (1 << bit);
  }

  // Flip the highest bit as well, based on the seed.
  out_bits ^= spacefillr::hash_u32(0, seed, 31) & (1 << 31);

  return out_bits;
}


//------------------------------------------------------

static inline uint32_t nested_uniform_scramble_base2_original_lk(uint32_t x, uint32_t seed) {
  x = reverse_bits(x);

  x += seed;
  x ^= x * 0x6c50b47cu;
  x ^= x * 0xb82f1e52u;
  x ^= x * 0xc7afe638u;
  x ^= x * 0x8d22f6e6u;

  x = reverse_bits(x);
  return x;
}

// void shuffled_scrambled_sobol4d_original_lk(uint32_t index, uint32_t seed,
//                                                    uint32_t X[4])
// {
//   index = nested_uniform_scramble_base2_original_lk(index, seed);
//   sobol4d(index, X);
//   for (int i = 0; i < 4; i++) {
//     X[i] = nested_uniform_scramble_base2_original_lk(X[i], hash_combine(seed, i));
//   }
// }

//------------------------------------------------------

static inline uint32_t nested_uniform_scramble_base2_v2(uint32_t x, uint32_t seed) {
  x = reverse_bits(x);

  x += seed;
  x ^= 0xdc967795;
  x *= 0x97b756bb;
  x ^= 0x866350b1;
  x *= 0x9e3779cd;

  x = reverse_bits(x);
  return x;
}

// void shuffled_scrambled_sobol4d_v2(uint32_t index, uint32_t seed,
//                                           uint32_t X[4])
// {
//   index = nested_uniform_scramble_base2_v2(index, seed);
//   sobol4d(index, X);
//   for (int i = 0; i < 4; i++) {
//     X[i] = nested_uniform_scramble_base2_v2(X[i], hash_combine(seed, i));
//   }
// }

//------------------------------------------------------

static inline uint32_t nested_uniform_scramble_base2_5round(uint32_t x, uint32_t seed) {
  x = reverse_bits(x);

  x *= 0x788aeeed;
  x ^= x * 0x41506a02;
  x += seed;
  x *= seed | 1;
  x ^= x * 0x7483dc64;

  x = reverse_bits(x);
  return x;
}

// void shuffled_scrambled_sobol4d_5round(uint32_t index, uint32_t seed,
//                                               uint32_t X[4])
// {
//   index = nested_uniform_scramble_base2_5round(index, seed);
//   sobol4d(index, X);
//   for (int i = 0; i < 4; i++) {
//     X[i] = nested_uniform_scramble_base2_5round(X[i], hash_combine(seed, i));
//   }
// }

//------------------------------------------------------

static inline uint32_t nested_uniform_scramble_base2_fast(uint32_t x, uint32_t seed) {
  x = reverse_bits(x);

  x += x << 2;
  x ^= x * 0xfe9b5742;
  x += seed;
  x *= seed | 1;

  x = reverse_bits(x);
  return x;
}

// void shuffled_scrambled_sobol4d_fast(uint32_t index, uint32_t seed,
//                                             uint32_t X[4])
// {
//   index = nested_uniform_scramble_base2_fast(index, seed);
//   sobol4d(index, X);
//   for (int i = 0; i < 4; i++) {
//     X[i] = nested_uniform_scramble_base2_fast(X[i], hash_combine(seed, i));
//   }
// }


namespace spacefillr {

static inline uint32_t sobol_u32(uint32_t index, uint32_t dimension, uint32_t scramble = 0) {
  if(dimension >= NumSobolDimensions) {
    throw std::runtime_error("Too many dimensions");
  }
  uint32_t v = scramble;
  for (int i = dimension * SobolMatrixSize; index != 0; index >>= 1, i++) {
    if (index & 1) {
      v ^= SobolMatrices32[i];
    }
  }
  return (v);
}


static inline  float u32_to_0_1_f32(uint32_t n) {
  return(std::fmin(n * 0x1p-32f /* 1/2^32 */,
                  FloatOneMinusEpsilon));
}


/// Scrambles `n` using fast hash-based Owen scrambling.
static inline  uint32_t owen_scramble_fast_u32(uint32_t x, uint32_t seed)  {
  x = reverse_bits(x);
  // Randomize the seed value.
  seed = hash_u32(seed, 0xa14a177d);

  // // Original Laine-Karras hash.
  // x = x.wrapping_add(seed);
  // x ^= x.wrapping_mul(0x6c50b47c);
  // x ^= x.wrapping_mul(0xb82f1e52);
  // x ^= x.wrapping_mul(0xc7afe638);
  // x ^= x.wrapping_mul(0x8d22f6e6);

  // // Fast, reasonable quality.
  // x = x + (x << 2);
  // x ^= x * (0xfe9b5742);
  // x = x + (seed);
  // x = x + (seed | 1);

  // // Medium-fast, best quality so far.
  x  = x * (0x788aeeed);
  x ^= x * (0x41506a02);
  x =  x + (seed);
  x =  x * (seed | 1);
  x ^= x * (0x7483dc64);

  return(reverse_bits(x));
}

/// Same as `sample()` except applies Owen scrambling using a fast hash-based
/// approach.
static inline float sobol_owen_fast_single(uint32_t index, uint32_t dimension, uint32_t seed) {
  return(u32_to_0_1_f32(owen_scramble_fast_u32(sobol_u32(index, dimension), seed)));
}

/// Same as `sample_owen_fast()` except it uses a slower "ground-truth"
/// implementation of Owen scrambling.
static inline float sobol_owen_reference_single(uint32_t index, uint32_t dimension, uint32_t seed)  {
  return(u32_to_0_1_f32(owen_scramble_u32(sobol_u32(index, dimension),seed)));
}

//No Scrambling
static inline float sobol_single(uint32_t index, uint32_t dimension, uint32_t scramble) {
  return(u32_to_0_1_f32(sobol_u32(index, dimension, scramble)));
}

//----------------------------------------------------------------------

} //namespace spacefillr

#endif
