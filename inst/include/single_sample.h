#ifndef SINGLESAMPLEH
#define SINGLESAMPLEH

#include "Rcpp.h"
#include "sobol.h"

namespace spacefillr {

static inline std::vector<float> sobol_fast_calc_std(uint64_t  i, unsigned int dim, unsigned int scramble) {
  std::vector<float> vals(dim);
  for(unsigned int j = 0; j < dim; j++) {
    vals[j] = spacefillr::sobol_owen_fast_single(i, j, scramble);
  }
  return(vals);
}

static inline Rcpp::List sobol_fast_calc_list(uint64_t  i, unsigned int dim, unsigned int scramble) {
  Rcpp::List vals(dim);
  for(unsigned int j = 0; j < dim; j++) {
    vals[j] = spacefillr::sobol_owen_fast_single(i, j, scramble);
  }
  return(vals);
}

}

#endif
