// Baseline.h: Pair of stations that together form a baseline (interferometer).
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// \file
/// \brief Pair of stations that together form a baseline (interferometer).

#ifndef DPPP_BASELINE_H
#define DPPP_BASELINE_H

#include <cstddef>
#include <utility>

namespace dp3 {
namespace base {

typedef std::pair<size_t, size_t> Baseline;

}  // namespace base
}  // namespace dp3

#endif
