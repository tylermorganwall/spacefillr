#ifndef RNGH
#define RNGH

#include <random>
#define extended extended_rng
#include "pcg/pcg_random.hpp"

class random_gen {
  public:
    random_gen(unsigned int seed) : rng(seed) {}
  random_gen() : rng(pcg_extras::seed_seq_from<std::random_device>{}) { }
  float unif_rand() {
    return(std::ldexp(rng(),-32));
  }
  uint32_t UniformUInt32(uint32_t b) {
    uint32_t threshold = (~b + 1u) % b;
    while (true) {
      uint32_t r = rng();
      if (r >= threshold)
        return r % b;
    }
  }
  pcg32 rng;

};


#endif
