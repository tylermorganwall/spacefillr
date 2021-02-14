/*
 * Copyright (C) Andrew Helmer 2020.
 * Licensed under MIT Open-Source License: see LICENSE.
 *
 * This file implements different methods of selecting the subquadrants in
 * between odd and even powers of 4 for the PMJ and PMJ02 algorithms. Compared
 * to random, they make a big difference for the overall error!
 */
#ifndef SAMPLE_GENERATION_SELECT_SUBQUAD_H_
#define SAMPLE_GENERATION_SELECT_SUBQUAD_H_

#include <iostream>
#include <utility>
#include <vector>

#include "util.h"

namespace pmj {
  typedef std::vector<std::pair<int, int>> (*subquad_fn)(
      const Point samples[], const int dim, random_gen& rng);

/*
* This will randomly choose once to swap X or swap Y, and will always swap
* X or Y for all subquadrants. For PMJ02, this ensures that the next set of
* samples are themselves a (0,2) sequence. Based off some basic analysis,
* it seems like this is the only way to maintain this property.
*
* Credit goes to Simon Brown for discovering this method with his Rust
* implementation: https://github.com/sjb3d/pmj
*/
std::vector<std::pair<int, int>> GetSubQuadrantsSwapXOrY(
 const Point samples[],
 const int dim, random_gen& rng);

/*
* Pick which subquadrants to use, using the ox-plowing technique from
* Christensen et al.
*/
std::vector<std::pair<int, int>> GetSubQuadrantsOxPlowing(
   const Point samples[],
   const int dim, random_gen& rng);

/*
* Pick which subquadrants to use randomly. No reason to actually use this:
* OxPlowing is better for pmj and ShuffleSwap is better for pmj(0,2).
*/
std::vector<std::pair<int, int>> GetSubQuadrantsRandomly(
   const Point samples[],
   const int dim, random_gen& rng);


std::vector<std::pair<int, int>> GetSubQuadrantsRandomly(
     const Point samples[],
                        const int dim, random_gen& rng) {
  const int quad_dim = dim / 2;
  const int n = quad_dim*quad_dim;

  std::vector<std::pair<int, int>> choices(n);

  for (int i = 0; i < n; i++) {
     const auto& sample = samples[i];
     int x_pos = sample.x * dim;
     int y_pos = sample.y * dim;

     if (UniformRand(0,1,rng) < 0.5) {
        choices[i] = {x_pos ^ 1, y_pos};
     } else {
        choices[i] = {x_pos, y_pos ^ 1};
     }
  }

  return choices;
}

std::vector<std::pair<int, int>> GetSubQuadrantsSwapXOrY(
     const Point samples[],
                        const int dim, random_gen& rng) {
  const int quad_dim = dim / 2;
  const int n = quad_dim*quad_dim;

  std::vector<std::pair<int, int>> choices(n);

  const bool swap_x = UniformRand(0,1,rng) < 0.5;

  for (int i = 0; i < n; i++) {
     const Point& sample = samples[i];
     int x_pos = sample.x * dim;
     int y_pos = sample.y * dim;

     if (swap_x) x_pos = x_pos ^ 1;
     else
        y_pos = y_pos ^ 1;

     choices[i] = {x_pos, y_pos};
  }

  return choices;
}

std::vector<std::pair<int, int>> GetSubQuadrantsOxPlowing(
     const Point samples[],
                        const int dim, random_gen& rng) {
  const int quad_dim = dim / 2;
  const int n = quad_dim*quad_dim;

  std::vector<std::pair<int, int>> choices(n);

  // First we want to get the subquadrant positions, and also the sampling order
  // from the original samples.
  std::vector<int> first_cells(n*2);
  std::vector<int> quadrant_order(n);
  for (int i = 0; i < n; i++) {
     const auto& sample = samples[i];
     int x_pos = sample.x * dim;
     int y_pos = sample.y * dim;
     const int quadrant_index = (y_pos / 2)*(quad_dim) + (x_pos / 2);
     first_cells[2*quadrant_index] = x_pos;
     first_cells[2*quadrant_index+1] = y_pos;
     quadrant_order[quadrant_index] = i;
  }

  // This method doesn't always work successfully, so we try a few times. In
  // the worst case, we'll always give a valid selection at the end anyway. In
  // practice, it virtually always succeeds within 10 attempts.
  for (int attempt = 0; attempt < 10; attempt++) {
     std::vector<int> choice_balance_x(quad_dim);
     std::vector<int> choice_balance_y(quad_dim);
     bool up = true;
     for (int col = 0; col < quad_dim; col++) {
        up = !up;
        for (int i = 0; i < quad_dim; i++) {
           const int row = up ? i : quad_dim - i - 1;

           const int quadrant_index = row*quad_dim + col;
           int x_pos = first_cells[2*quadrant_index];
           int y_pos = first_cells[2*quadrant_index+1];

           const bool last = (i == quad_dim - 1);
           const int balance_y = choice_balance_y[row];
           const int balance_x = choice_balance_x[col];

           bool swap_x = false;
           if (balance_y != 0 && !last) {
              swap_x = (balance_y > 0) != (y_pos & 1);
           } else if (balance_x != 0) {
              swap_x = (balance_x > 0) == (x_pos & 1);
           } else {
              swap_x = UniformRand(0,1,rng) < 0.5;
           }

           x_pos = swap_x ? x_pos ^ 1 : x_pos;
           y_pos = (!swap_x) ? y_pos ^ 1 : y_pos;

           choices[quadrant_order[quadrant_index]].first = x_pos;
           choices[quadrant_order[quadrant_index]].second = y_pos;

           choice_balance_x[col] += (x_pos & 1) ? 1 : -1;
           choice_balance_y[row] += (y_pos & 1) ? 1 : -1;
        }
     }

     // Always unbalanced, just return.
     if (n == 1) {
        return choices;
     }

     bool is_balanced = true;
     for (int row = 0; row < quad_dim; row++) {
        if (choice_balance_y[row] != 0) {
           is_balanced = false;
           break;
        }
     }
     if (is_balanced) break;
  }

  return choices;
}

}  // namespace pmj

#endif  // SAMPLE_GENERATION_SELECT_SUBQUAD_H_
