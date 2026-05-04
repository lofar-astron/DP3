// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_BASE_CALTYPE_H_
#define DP3_BASE_CALTYPE_H_

#include <string>

namespace dp3::base {

enum class CalType {
  kScalar,
  kScalarAmplitude,
  kScalarPhase,
  kDiagonal,
  kDiagonalAmplitude,
  kDiagonalPhase,
  kFullJones,
  kTec,
  kTecAndDelay,
  kTecAndPhase,
  kTecPhaseAndDelay,
  kTecScreen,
  kRotationAndDiagonal,
  kRotation,
  kFaradayRotation,
  kLeakage,
  kLeakageAmplitude
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
    case CalType::kLeakage:
    case CalType::kLeakageAmplitude:
      return 4;
    case CalType::kScalar:
    case CalType::kScalarAmplitude:
    case CalType::kScalarPhase:
    case CalType::kTec:
    case CalType::kTecAndDelay:
    case CalType::kTecAndPhase:
    case CalType::kTecPhaseAndDelay:
    case CalType::kTecScreen:
      return 1;
  }
  return 0;
}

}  // namespace dp3::base

#endif  // DP3_CALTYPE_H
