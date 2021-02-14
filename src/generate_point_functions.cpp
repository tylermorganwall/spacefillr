#include <Rcpp.h>
using namespace Rcpp;

#include "pj.h"
#include "pmj.h"

#include "pmj02.h"
#include "sobol.h"
#include "halton_sampler.h"

#include "rng.h"

#include "util.h"

// [[Rcpp::export]]
List rcpp_generate_sobol_set(unsigned long long  N, unsigned int dim, unsigned int scramble) {
  List final_set(N*dim);
  int counter = 0;
  for(unsigned int j = 0; j < dim; j++) {
    for(unsigned long long i = 0; i < N; i++) {
      final_set(counter) = spacefillr::sobol_single(i, j, scramble);
      counter = counter + 1;
    }
  }
  return(final_set);
}

// [[Rcpp::export]]
List rcpp_generate_sobol_owen_set(unsigned long long  N, unsigned int dim, unsigned int scramble) {
  List final_set(N*dim);
  int counter = 0;
  for(unsigned int j = 0; j < dim; j++) {
    for(unsigned long long i = 0; i < N; i++) {
      final_set(counter) = spacefillr::sobol_owen_reference_single(i, j, scramble);
      counter = counter + 1;
    }
  }
  return(final_set);
}

// [[Rcpp::export]]
List rcpp_generate_sobol_owen_fast_set(unsigned long long  N, unsigned int dim, unsigned int scramble) {
  List final_set(N*dim);
  int counter = 0;
  for(unsigned int j = 0; j < dim; j++) {
    for(unsigned long long i = 0; i < N; i++) {
      final_set(counter) = spacefillr::sobol_owen_fast_single(i, j, scramble);
      counter = counter + 1;
    }
  }
  return(final_set);
}

// [[Rcpp::export]]
double rcpp_generate_sobol_owen_fast_single(unsigned long long  i, unsigned int dim, unsigned int scramble) {
  return(spacefillr::sobol_owen_fast_single(i, dim, scramble));
}


// [[Rcpp::export]]
List rcpp_generate_halton_faure_set(unsigned long long  N, unsigned int dim) {
  List final_set(N*dim);
  Halton_sampler hs;
  hs.init_faure();
  int counter = 0;
  for(unsigned int j = 0; j < dim; j++) {
    for(unsigned long long i = 0; i < N; i++) {
      final_set(counter) = hs.sample(j,i);
      counter = counter + 1;
    }
  }
  return(final_set);
}



// [[Rcpp::export]]
List rcpp_generate_halton_random_set(unsigned long long  N, unsigned int dim, unsigned int seed) {
  List final_set(N*dim);
  random_gen rng(seed);
  Halton_sampler hs;
  hs.init_random(rng.rng);
  int counter = 0;
  for(unsigned int j = 0; j < dim; j++) {
    for(unsigned long long i = 0; i < N; i++) {
      final_set(counter) = hs.sample(j,i);
      counter = counter + 1;
    }
  }
  return(final_set);
}


// [[Rcpp::export]]
double rcpp_generate_halton_faure_single(unsigned long long  i, unsigned int dim) {
  Halton_sampler hs;
  hs.init_faure();
  return(hs.sample(i,dim));
}



// [[Rcpp::export]]
double rcpp_generate_halton_random_single(unsigned long long  i, unsigned int dim, unsigned int seed) {
  random_gen rng(seed);
  Halton_sampler hs;
  hs.init_random(rng.rng);
  return(hs.sample(i,dim));
}

//PJ + PMJ

// [[Rcpp::export]]
List rcpp_generate_pj_set(unsigned long long  N, int seed) {
  List final_set(N*2);
  random_gen rng(seed);
  std::unique_ptr<pmj::Point[]> points = pmj::GetProgJitteredSamples(N, rng);
  int counter = 0;
  for(int i = 0; i < N; i++) {
    final_set(counter) = points[i].x;
    final_set(counter+1) = points[i].y;
    counter += 2;
  }
  return(final_set);
}

// [[Rcpp::export]]
List rcpp_generate_pmj_set(unsigned long long  N, int seed) {
  List final_set(N*2);
  random_gen rng(seed);
  std::unique_ptr<pmj::Point[]> points = pmj::GetProgMultiJitteredSamples(N, rng);
  int counter = 0;
  for(int i = 0; i < N; i++) {
    final_set(counter) = points[i].x;
    final_set(counter+1) = points[i].y;
    counter += 2;
  }
  return(final_set);
}


// [[Rcpp::export]]
List rcpp_generate_pmjbn_set(unsigned long long  N, int seed) {
  List final_set(N*2);
  random_gen rng(seed);
  std::unique_ptr<pmj::Point[]> points = pmj::GetProgMultiJitteredSamplesWithBlueNoise(N, rng);
  int counter = 0;
  for(int i = 0; i < N; i++) {
    final_set(counter) = points[i].x;
    final_set(counter+1) = points[i].y;
    counter += 2;
  }
  return(final_set);
}

// [[Rcpp::export]]
List rcpp_generate_pmj02_set(unsigned long long  N, int seed) {
  List final_set(N*2);
  random_gen rng(seed);
  std::unique_ptr<pmj::Point[]> points = pmj::GetPMJ02Samples(N, rng);
  int counter = 0;
  for(int i = 0; i < N; i++) {
    final_set(counter) = points[i].x;
    final_set(counter+1) = points[i].y;
    counter += 2;
  }
  return(final_set);
}

// [[Rcpp::export]]
List rcpp_generate_pmj02bn_set(unsigned long long  N, int seed) {
  List final_set(N*2);
  random_gen rng(seed);
  std::unique_ptr<pmj::Point[]> points = pmj::GetPMJ02SamplesWithBlueNoise(N, rng);
  int counter = 0;
  for(int i = 0; i < N; i++) {
    final_set(counter) = points[i].x;
    final_set(counter+1) = points[i].y;
    counter += 2;
  }
  return(final_set);
}





