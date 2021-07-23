// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_CALTYPE_H
#define DP3_CALTYPE_H

#include <string>

namespace dp3 {
namespace base {

enum class CalType {
  kScalar,
  kScalarAmplitude,
  kScalarPhase,
  kDiagonal,
  kDiagonalAmplitude,
  kDiagonalPhase,
  kFullJones,
  kTecAndPhase,
  kTec,
  kTecScreen,
  kRotationAndDiagonal,
  kRotation
};

/// Convert string to a CalType
CalType StringToCalType(const std::string& mode);

/// Convert CalType to a string
std::string ToString(CalType caltype);

}  // namespace base
}  // namespace dp3

#endif  // DP3_CALTYPE_H