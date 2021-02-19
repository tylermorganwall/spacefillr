/*
 * Copyright (C) Andrew Helmer 2020.
 * Licensed under MIT Open-Source License: see LICENSE.
 *
 * Generate Progressive Multi-Jittered Sequences from
 * "Progressive Multi-Jittered Sample Sequences", Christensen et al. 2018. The
 * non-best-candidate sequences generate about 2 million samples/sec for me,
 * it's not very optimized compared to Christensen's paper.
 * The best candidate sequences only do about ~230k samples/sec.
 *
 * If you're reading this code for the first time and want to understand the
 * algorithm, start with the function "GenerateSamples".
 *
 */
#ifndef SAMPLE_GENERATION_PMJ_H_
#define SAMPLE_GENERATION_PMJ_H_

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <random>
#include <utility>
#include <vector>

#include "util.h"
#include "select_subquad.h"

namespace pmj {


// Generates progressive multi-jittered samples without blue noise properties.
// Takes in a number of samples.
std::unique_ptr<Point[]> GetProgMultiJitteredSamples(
    const int num_samples, random_gen& rng);

// Generates progressive multi-jittered samples with blue noise properties, i.e.
// using best-candidate points.
std::unique_ptr<Point[]> GetProgMultiJitteredSamplesWithBlueNoise(
    const int num_samples, random_gen& rng);

/*
 * The SampleSet is a class that contains the generated samples, as well as the
 * currently populated strata. It's used to generate new samples within the
 * unpopulated strata.
 */
class SampleSet2 {
    public:
        explicit SampleSet2(const int num_samples,
                           const int num_candidates,
                           random_gen& rng_)
            : num_candidates_(num_candidates), rng(rng_) {
            samples_ = std::unique_ptr<Point[]>(new Point[num_samples]());

            int grid_memory_size = 1;
            while (grid_memory_size < num_samples)
                grid_memory_size <<= 2;
            sample_grid_ = std::unique_ptr<const Point*[]>(new const Point*[grid_memory_size]());
            x_strata_.resize(grid_memory_size);
            y_strata_.resize(grid_memory_size);
        }

        // This generates a new sample at the current index, given the X position
        // and Y position of the subquadrant. It won't generate a new sample in an
        // existing strata.
        void GenerateNewSample(const int sample_index,
                               const int x_pos,
                               const int y_pos);

        // This function should be called after every power of 2 samples.
        void SubdivideStrata();

        // Get all the samples at the end.
        std::unique_ptr<Point[]> ReleaseSamples() {
            return std::move(samples_);
        }

        const Point& sample(const int sample_index) const {
            return samples_[sample_index];
        }
        const Point* samples() const {
            return samples_.get();
        }
        const int dim() const { return dim_; }

        private:
            // Generates a valid sample at the given cell position, subject to
            // stratification.
            Point GetCandidateSample(const int x_pos, const int y_pos);

            // Adds a new point at index i. Updates the necessary data structures.
            void AddSample(const int i, const Point& sample);

            std::unique_ptr<Point[]> samples_;

            // This could be SIGNIFICANTLY optimized, especially for best-candidate
            // sampling, by actually storing these as a binary tree, rather than a linear
            // array, where each node represents whether there are any unoccupied strata
            // in a range. Then for each sample, you could traverse the tree to get the
            // unoccupied strata. This would be similar to Matt Pharr's optimization for
            // PMJ02.
            std::vector<bool> x_strata_ {false};
            std::vector<bool> y_strata_ {false};

            // The sample grid is used for nearest neighbor lookups.
            std::unique_ptr<const Point*[]> sample_grid_;

            int n_ = 1;  // Number of samples in the next pass.
            bool is_power_of_4_ = true;  // Whether n is a power of 4.
            int dim_ = 1;  // Number of cells in one dimension in next pass, i.e. sqrt(n).
            double grid_size_ = 1.0;  // 1.0 / dim_

            // Number of candidates to use for best-candidate sampling.
            const int num_candidates_;
            random_gen rng;
};

// This generates a sample within the grid position, verifying that it doesn't
// overlap strata with any other sample.
double Get1DStrataSample(const int pos,
                         const int n,
                         const double grid_size,
                         const std::vector<bool>& strata,
                         random_gen& rng) {
    while (true) {
        double val = UniformRand(pos*grid_size, (pos+1)*grid_size, rng);
        int strata_pos = val * n;
        if (!strata[strata_pos]) {
            return val;
        }
    }
}

Point SampleSet2::GetCandidateSample(const int x_pos,
                                    const int y_pos) {
    return {Get1DStrataSample(x_pos, n_, grid_size_, x_strata_, rng),
            Get1DStrataSample(y_pos, n_, grid_size_, y_strata_, rng)};
}

void SampleSet2::GenerateNewSample(const int sample_index,
                                  const int x_pos,
                                  const int y_pos) {
    Point best_candidate;
    if (num_candidates_ <= 1) {
        best_candidate = GetCandidateSample(x_pos, y_pos);
    } else {
        std::vector<Point> candidate_samples(num_candidates_);
        for (int i = 0; i < num_candidates_; i++) {
            candidate_samples[i] = GetCandidateSample(x_pos, y_pos);
        }

        best_candidate = GetBestCandidateOfSamples(
            candidate_samples, sample_grid_.get(), dim_);
    }
    AddSample(sample_index, best_candidate);
}

void SampleSet2::SubdivideStrata() {
    const int old_n = n_;

    n_ *= 2;
    is_power_of_4_ = !is_power_of_4_;
    if (!is_power_of_4_) {
        dim_ *= 2;
        grid_size_ *= 0.5;
    }

    std::fill_n(sample_grid_.get(), n_, nullptr);
    std::fill_n(x_strata_.begin(), n_, 0);
    std::fill_n(y_strata_.begin(), n_, 0);
    for (int i = 0; i < old_n; i++) {
        const auto& sample = samples_[i];

        x_strata_[sample.x * n_] = true;
        y_strata_[sample.y * n_] = true;

        const int x_pos = sample.x * dim_, y_pos = sample.y * dim_;
        sample_grid_[y_pos*dim_ + x_pos] = &sample;
    }
}

void SampleSet2::AddSample(const int i,
                          const Point& sample) {
    samples_[i] = sample;

    x_strata_[sample.x * n_] = true;
    y_strata_[sample.y * n_] = true;

    const int x_pos = sample.x * dim_, y_pos = sample.y * dim_;
    sample_grid_[y_pos*dim_ + x_pos] = &(samples_[i]);
}

/*
 * The core of Christensen et al.'s algorithm.
 */
std::unique_ptr<Point[]> GenerateSamples2(
        const int num_samples,
        const int num_candidates,
        random_gen& rng,
        const subquad_fn subquad_func = &GetSubQuadrantsOxPlowing) {
    SampleSet2 sample_set(num_samples, num_candidates, rng);

    // Generate first sample.
    sample_set.GenerateNewSample(0, 0, 0);

    int quadrants = 1;
    while (quadrants < num_samples) {
        sample_set.SubdivideStrata();

        // For every sample, we first generate the diagonally opposite one at the
        // current grid level.
        for (int i = 0;
             i < quadrants && quadrants+i < num_samples;
             i++) {
            const auto& sample = sample_set.sample(i);

            int x_pos = sample.x * sample_set.dim();
            int y_pos = sample.y * sample_set.dim();

            sample_set.GenerateNewSample(quadrants+i, x_pos ^ 1, y_pos ^ 1);
            if (quadrants+i >= num_samples) {
                break;
            }
        }

        if (2*quadrants >= num_samples) break;

        // Now we generate samples in the remaining subquadrants.
        sample_set.SubdivideStrata();

        // We want to make balanced choices here regarding which subquadrants to
        // use, so we precompute them in a special way.
        auto sub_quad_choices =
            (*subquad_func)(sample_set.samples(), sample_set.dim(), rng);
        // Iterate over the quadrants and generate a new point in each quadrant.
        for (int i = 0;
             i < quadrants && 2*quadrants+i < num_samples;
             i++) {
            sample_set.GenerateNewSample(
                2*quadrants+i, sub_quad_choices[i].first, sub_quad_choices[i].second);
        }

        for (int i = 0;
             i < quadrants && 3*quadrants+i < num_samples;
             i++) {
            // Get the one diagonally opposite to the one we just got.
            sample_set.GenerateNewSample(3*quadrants+i,
                                         sub_quad_choices[i].first ^ 1,
                                         sub_quad_choices[i].second ^ 1);
        }

        quadrants *= 4;
    }

    return sample_set.ReleaseSamples();
}

std::unique_ptr<Point[]> GetProgMultiJitteredSamples(
        const int num_samples, random_gen& rng) {
    return GenerateSamples2(num_samples, /*num_candidates=*/1, rng);
}
std::unique_ptr<Point[]> GetProgMultiJitteredSamplesWithBlueNoise(
        const int num_samples, random_gen& rng) {
    return GenerateSamples2(num_samples, kBestCandidateSamples, rng);
}

}  // namespace pmj

#endif  // SAMPLE_GENERATION_PMJ_H_
