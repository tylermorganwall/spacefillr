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

#include <array>
#include <memory>
#include <utility>
#include <vector>

#include "util.h"

namespace pmj {

// Progressive jittered samples shouldn't really be used, it's more just a
// learning example.
std::unique_ptr<pmj::Point[]> GetProgJitteredSamples(
    const int num_samples, random_gen& rng);

}  // namespace pmj

#endif  // SAMPLE_GENERATION_PJ_H_
