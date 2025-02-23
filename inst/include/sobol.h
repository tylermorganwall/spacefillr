#ifndef SOBOLH
#define SOBOLH

#include "sobolmatrices.h"
#include "sobol_directions.h"
#include "siphash.h"
#include <stdexcept>
#include <cmath>

static const float FloatOneMinusEpsilon = 0.99999994;

static inline  float u32_to_0_1_f32(uint32_t n) {
  return(std::fmin(n * 0x1p-32f /* 1/2^32 */,
                   FloatOneMinusEpsilon));
}
static inline uint32_t hash_combine(uint32_t seed, uint32_t v) {
  return seed ^ (v + (seed << 6) + (seed >> 2));
}

static inline uint32_t hash(uint32_t x) {
  // finalizer from murmurhash3
  x ^= x >> 16;
  x *= 0x85ebca6bu;
  x ^= x >> 13;
  x *= 0xc2b2ae35u;
  x ^= x >> 16;
  return x;
}

static inline uint32_t sobol(uint32_t index, uint32_t dim) {
  uint32_t X = 0;
  for (int bit = 0; bit < 32; bit++) {
    int mask = (index >> bit) & 1;
    X ^= mask * SPACEFILLR_SOBOL_DIRECTIONS[dim][bit];
  }
  return X;
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


namespace spacefillr {


/// Scrambles `n` using fast hash-based Owen scrambling.
static inline  uint32_t owen_scramble_fast_u32(uint32_t x, uint32_t seed)  {
  x = reverse_bits(x);
  // Randomize the seed value.
  seed = hash_u32(seed, 0xa14a177d);

  x ^= x * 0x3d20adea;
  x += seed;
  x *= (seed >> 16) | 1;
  x ^= x * 0x05526c56;
  x ^= x * 0x53a22864;

  return(reverse_bits(x));
}


static inline uint32_t sobol_u32(uint32_t index, uint32_t dimension, uint32_t scramble = 0) {
  if(dimension >= NumSobolDimensions) {
    throw std::runtime_error("Too many dimensions");
  }
  index = owen_scramble_fast_u32(index, scramble);
  uint32_t v = 0;
  for (int i = dimension * SobolMatrixSize; index != 0; index >>= 1, i++) {
    if (index & 1) {
      v ^= SobolMatrices32[i];
    }
  }
  return (v);
}

/// Same as `sample()` except applies Owen scrambling using a fast hash-based
/// approach.
static inline float sobol_owen_single(uint32_t index, uint32_t dimension, uint32_t seed) {
  if(dimension > 21201) {
    throw std::runtime_error("Too many dimensions");
  }
  // uint32_t i = owen_scramble_fast_u32(index, seed);
  return(u32_to_0_1_f32(owen_scramble_fast_u32(sobol(owen_scramble_fast_u32(index, seed), dimension),
                                hash_combine(seed, dimension))));
}

//No Scrambling
static inline float sobol_single(uint32_t index, uint32_t dimension, uint32_t scramble) {
  return(u32_to_0_1_f32(sobol_u32(index, dimension, scramble)));
}

//----------------------------------------------------------------------

} //namespace spacefillr

#endif
