// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_BASE_CALTYPE_H_
#define DP3_BASE_CALTYPE_H_

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
  kRotation,
  kFaradayRotation
};

/// Convert string to a CalType
CalType StringToCalType(const std::string& mode);

/// Convert CalType to a string
std::string ToString(CalType caltype);

constexpr size_t GetNPolarizations(CalType cal_type) {
  switch (cal_type) {
    case CalType::kDiagonal:
    case CalType::kDiagonalPhase:
    case CalType::kDiagonalAmplitude:
      return 2;
    case CalType::kFullJones:
    case CalType::kRotationAndDiagonal:
    case CalType::kRotation:
    case CalType::kFaradayRotation:
      return 4;
    case CalType::kScalar:
    case CalType::kScalarAmplitude:
    case CalType::kScalarPhase:
    case CalType::kTecAndPhase:
    case CalType::kTec:
    case CalType::kTecScreen:
      return 1;
  }
  return 0;
}

}  // namespace base
}  // namespace dp3

#endif  // DP3_CALTYPE_H
