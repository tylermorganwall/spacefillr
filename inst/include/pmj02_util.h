/*
 * Copyright (C) Andrew Helmer 2020.
 * Licensed under MIT Open-Source License: see LICENSE.
 *
 * Really just implements Matt Pharr's algorithm for speeding up PMJ(0,2).
 */
#ifndef SAMPLE_GENERATION_PMJ02_UTIL_H_
#define SAMPLE_GENERATION_PMJ02_UTIL_H_

#include <utility>
#include <vector>

namespace pmj {

using std::vector;

// Implementation of "Efficient Generation of Points that Satisfy
// Two-Dimensional Elementary Intervals" by Matt Pharr (2019). Given a set of
// strata corresponding to all elementary (0,2) intervals, and a set of
// positions on a square grid, returns the narrowest X strata and Y strata that
// are unoccupied.
std::pair<std::vector<int>, std::vector<int>>
    GetValidStrata(const int x_pos,
                   const int y_pos,
                   const std::vector<std::vector<bool>>& strata);


void GetXStrata(const int x_pos,
                const int y_pos,
                const int strata_index,
                const vector<vector<bool>>& strata,
                vector<int>* x_strata) {
    const int strata_n_cols = 1 << (strata.size() - strata_index - 1);
    const bool is_occupied =
        strata[strata_index][y_pos*strata_n_cols + x_pos];

    if (!is_occupied) {
        if (strata_index == 0) {
            // We're at the Nx1 leaf.
            x_strata->push_back(x_pos);
        } else {
            GetXStrata(x_pos * 2, y_pos / 2, strata_index - 1, strata, x_strata);
            GetXStrata(x_pos * 2 + 1, y_pos / 2, strata_index - 1, strata, x_strata);
        }
    }
}
void GetYStrata(const int x_pos,
                const int y_pos,
                const int strata_index,
                const vector<vector<bool>>& strata,
                vector<int>* y_strata) {
    const int strata_n_cols = 1 << (strata.size() - strata_index - 1);
    const bool is_occupied =
        strata[strata_index][y_pos*strata_n_cols + x_pos];

    if (!is_occupied) {
        if (strata_n_cols == 1) {
            // We're at the 1xN leaf.
            y_strata->push_back(y_pos);
        } else {
            GetYStrata(x_pos / 2, y_pos * 2, strata_index + 1, strata, y_strata);
            GetYStrata(x_pos / 2, y_pos * 2 + 1, strata_index + 1, strata, y_strata);
        }
    }
}

std::pair<vector<int>, vector<int>> GetValidStrata(
        const int x_pos, const int y_pos, const vector<vector<bool>>& strata) {
    std::pair<vector<int>, vector<int>> valid_strata = {{}, {}};

    if (strata.size() % 2 == 1) {
        GetXStrata(x_pos, y_pos, strata.size()/2, strata, &valid_strata.first);
        GetYStrata(x_pos, y_pos, strata.size()/2, strata, &valid_strata.second);
    } else {
        GetXStrata(x_pos, y_pos/2, strata.size()/2-1, strata, &valid_strata.first);
        GetYStrata(x_pos/2, y_pos, strata.size()/2, strata, &valid_strata.second);
    }

    return valid_strata;
}

}  // namespace pmj

#endif  // SAMPLE_GENERATION_PMJ02_UTIL_H_
