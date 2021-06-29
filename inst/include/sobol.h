#ifndef SOBOLH
#define SOBOLH

#include "sobolmatrices.h"
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


static uint32_t directions[5][32] = {
  {0x80000000, 0x40000000, 0x20000000, 0x10000000,
  0x08000000, 0x04000000, 0x02000000, 0x01000000,
  0x00800000, 0x00400000, 0x00200000, 0x00100000,
  0x00080000, 0x00040000, 0x00020000, 0x00010000,
  0x00008000, 0x00004000, 0x00002000, 0x00001000,
  0x00000800, 0x00000400, 0x00000200, 0x00000100,
  0x00000080, 0x00000040, 0x00000020, 0x00000010,
  0x00000008, 0x00000004, 0x00000002, 0x00000001},

  {0x80000000, 0xc0000000, 0xa0000000, 0xf0000000,
  0x88000000, 0xcc000000, 0xaa000000, 0xff000000,
  0x80800000, 0xc0c00000, 0xa0a00000, 0xf0f00000,
  0x88880000, 0xcccc0000, 0xaaaa0000, 0xffff0000,
  0x80008000, 0xc000c000, 0xa000a000, 0xf000f000,
  0x88008800, 0xcc00cc00, 0xaa00aa00, 0xff00ff00,
  0x80808080, 0xc0c0c0c0, 0xa0a0a0a0, 0xf0f0f0f0,
  0x88888888, 0xcccccccc, 0xaaaaaaaa, 0xffffffff},

  {0x80000000, 0xc0000000, 0x60000000, 0x90000000,
  0xe8000000, 0x5c000000, 0x8e000000, 0xc5000000,
  0x68800000, 0x9cc00000, 0xee600000, 0x55900000,
  0x80680000, 0xc09c0000, 0x60ee0000, 0x90550000,
  0xe8808000, 0x5cc0c000, 0x8e606000, 0xc5909000,
  0x6868e800, 0x9c9c5c00, 0xeeee8e00, 0x5555c500,
  0x8000e880, 0xc0005cc0, 0x60008e60, 0x9000c590,
  0xe8006868, 0x5c009c9c, 0x8e00eeee, 0xc5005555},

  {0x80000000, 0xc0000000, 0x20000000, 0x50000000,
  0xf8000000, 0x74000000, 0xa2000000, 0x93000000,
  0xd8800000, 0x25400000, 0x59e00000, 0xe6d00000,
  0x78080000, 0xb40c0000, 0x82020000, 0xc3050000,
  0x208f8000, 0x51474000, 0xfbea2000, 0x75d93000,
  0xa0858800, 0x914e5400, 0xdbe79e00, 0x25db6d00,
  0x58800080, 0xe54000c0, 0x79e00020, 0xb6d00050,
  0x800800f8, 0xc00c0074, 0x200200a2, 0x50050093},

  {0x80000000, 0x40000000, 0x20000000, 0xb0000000,
  0xf8000000, 0xdc000000, 0x7a000000, 0x9d000000,
  0x5a800000, 0x2fc00000, 0xa1600000, 0xf0b00000,
  0xda880000, 0x6fc40000, 0x81620000, 0x40bb0000,
  0x22878000, 0xb3c9c000, 0xfb65a000, 0xddb2d000,
  0x78022800, 0x9c0b3c00, 0x5a0fb600, 0x2d0ddb00,
  0xa2878080, 0xf3c9c040, 0xdb65a020, 0x6db2d0b0,
  0x800228f8, 0x400b3cdc, 0x200fb67a, 0xb00ddb9d}
};

static inline uint32_t sobol(uint32_t index, uint32_t dim) {
  if (dim > 4) return 0;
  uint32_t X = 0;
  for (int bit = 0; bit < 32; bit++) {
    int mask = (index >> bit) & 1;
    X ^= mask * directions[dim][bit];
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

  x ^= x * 0x3d20adea;
  x += seed;
  x *= (seed >> 16) | 1;
  x ^= x * 0x05526c56;
  x ^= x * 0x53a22864;

  // // Medium-fast, best quality so far.
  // x  = x * (0x788aeeed);
  // x ^= x * (0x41506a02);
  // x =  x + (seed);
  // x =  x * (seed | 1);
  // x ^= x * (0x7483dc64);

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
