/*
 * Copyright (C) Andrew Helmer 2020.
 * Licensed under MIT Open-Source License: see LICENSE.
 *
 * Implements the most basic algorithm from Christensen et al., Progressive
 * Jittered Sampling. There's really no reason to use this, but the code could
 * be instructional.
 */
#ifndef SAMPLE_GENERATION_PJ_H_
#define SAMPLE_GENERATION_PJ_H_

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <random>
#include <utility>
#include <vector>

#include "util.h"
#include "rng.h"

namespace pmj {

// Progressive jittered samples shouldn't really be used, it's more just a
// learning example.
std::unique_ptr<pmj::Point[]> GetProgJitteredSamples(
    const int num_samples, random_gen& rng);

Point RandomSample(
        double min_x, double max_x, double min_y, double max_y, random_gen& rng) {
    return {UniformRand(min_x, max_x, rng), UniformRand(min_y, max_y, rng)};
}

Point GetSample(
        const int x_pos, const int y_pos, const double grid_size, random_gen& rng) {
    return RandomSample(x_pos*grid_size, (x_pos+1)*grid_size,
                        y_pos*grid_size, (y_pos+1)*grid_size, rng);
}

void GenerateSamplesForQuadrant(
        const Point& sample,
        const int num_samples,
        const int n,
        const int i,
        const int x_pos,
        const int y_pos,
        const double grid_size,
        Point* samples,
        random_gen& rng) {
    // Generate diagonally opposite.
    samples[n+i] = GetSample(x_pos ^ 1, y_pos ^ 1, grid_size, rng);

    if (2*n+i >= num_samples) {
        return;
    }

    // Pick one of the two adjacent cells to generate new sample.
    int new_x_pos = x_pos;
    int new_y_pos = y_pos;
    if (UniformRand(0, 1, rng) < 0.5) {
        new_x_pos = x_pos ^ 1;
    } else {
        new_y_pos = y_pos ^ 1;
    }
    samples[2*n+i] = GetSample(new_x_pos, new_y_pos, grid_size, rng);

    // Generate a sample in the diagonal of the previous cell.
    if (3*n+i >= num_samples) {
        return;
    }

    samples[3*n+i] = GetSample(new_x_pos ^ 1, new_y_pos ^ 1, grid_size, rng);
}

std::unique_ptr<Point[]> GetProgJitteredSamples(
        const int num_samples, random_gen& rng) {
    auto samples = std::unique_ptr<Point[]>(new Point[num_samples]());

    // Generate first sample randomly.
    samples[0] = RandomSample(0, 1, 0, 1, rng);

    int n = 1;  // Number of samples in previous pass.
    int dim = 2;  // The number of subquadrants in one dimension.
    double grid_size = 0.5;  // The subquadrant size in one dimension, 1.0 / dim.
    while (n < num_samples) {
        for (int i = 0; i < n && n+i < num_samples; i++) {
            const auto& sample = samples[i];

            int x_pos = sample.x * dim;
            int y_pos = sample.y * dim;

            GenerateSamplesForQuadrant(
                sample, num_samples, n, i, x_pos, y_pos, grid_size, samples.get(), rng);
        }
        n *= 4;
        dim *= 2;
        grid_size *= 0.5;
    }

    return samples;
}

} //namespace pmj

#endif  // SAMPLE_GENERATION_PJ_H_
