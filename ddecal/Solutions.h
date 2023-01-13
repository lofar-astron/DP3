// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_DDECAL_SOLUTIONS_H_
#define DP3_DDECAL_SOLUTIONS_H_

#include <aocommon/xt/span.h>
#include <xtensor/xtensor.hpp>

namespace dp3 {
namespace ddecal {
using SolutionsTensor = xt::xtensor<std::complex<double>, 4>;
using SolutionsSpan = aocommon::xt::Span<std::complex<double>, 4>;
}  // namespace ddecal
}  // namespace dp3

#endif  // DP3_DDECAL_SOLUTIONS_H_