#ifndef RNGH
#define RNGH

#include <random>
#define extended extended_rng
#include "pcg/pcg_random.hpp"

class random_gen {
public:
  // Define result_type for compatibility with standard library requirements
  using result_type = uint32_t;

  random_gen(unsigned int seed) : rng(seed) {}
  random_gen() : rng(pcg_extras::seed_seq_from<std::random_device>{}) {}

  // The uniform random number generator function
  float unif_rand() {
    return std::ldexp(rng(), -32);
  }

  uint32_t UniformUInt32(uint32_t b) {
    uint32_t threshold = (~b + 1u) % b;
    while (true) {
      uint32_t r = rng();
      if (r >= threshold)
        return r % b;
    }
  }

  // Call operator to satisfy UniformRandomBitGenerator requirements
  result_type operator()() {
    return rng();
  }

  // Static methods to satisfy UniformRandomBitGenerator requirements
  static constexpr result_type min() {
    return pcg32::min();
  }

  static constexpr result_type max() {
    return pcg32::max();
  }

  pcg32 rng;
};

#endif
