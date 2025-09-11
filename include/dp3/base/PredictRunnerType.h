// PredictRunnerType.h: Header for selecting between OnePredict and FastPredict
// Copyright (C) 2025 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_BASE_PREDICTRUNNERTYPE_H_
#define DP3_BASE_PREDICTRUNNERTYPE_H_

#include <steps/OnePredict.h>
#include <steps/FastPredict.h>

namespace dp3::base {

// Use FastPredict if USE_FAST_PREDICT is defined in the user's CMake
// configuration.
#ifdef USE_FAST_PREDICT
using PredictRunnerType = steps::FastPredict;
#else
using PredictRunnerType = steps::OnePredict;
#endif

}  // namespace dp3::base

#endif  // DP3_BASE_PREDICTRUNNERTYPE_H_