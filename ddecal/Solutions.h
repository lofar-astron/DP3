// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_DDECAL_SOLUTIONS_H_
#define DP3_DDECAL_SOLUTIONS_H_

#include <aocommon/xt/span.h>
#include <xtensor/xtensor.hpp>

namespace dp3 {
namespace ddecal {
/// Data structures for solutions. The dimensions are
/// (nr. channel blocks, nr. antennas, nr. solutions, nr. polarizations).
/// The number of solutions (third dimension) depends on the number of
/// directions and the number of solutions per direction.
/// Different directions may have different solution counts.
/// @sa dp3::ddecal::SolveData::ChannelBlockData::NSolutionsForDirection
/// @{
using SolutionTensor = xt::xtensor<std::complex<double>, 4>;
using SolutionSpan = aocommon::xt::Span<std::complex<double>, 4>;
/// @}
}  // namespace ddecal
}  // namespace dp3

#endif  // DP3_DDECAL_SOLUTIONS_H_