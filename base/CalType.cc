// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "CalType.h"

#include <stdexcept>

namespace dp3 {
namespace base {

CalType StringToCalType(const std::string& modestr) {
  if (modestr == "diagonal" || modestr == "complexgain")
    return CalType::kDiagonal;
  else if (modestr == "diagonalphase" || modestr == "phaseonly")
    return CalType::kDiagonalPhase;
  else if (modestr == "diagonalamplitude" || modestr == "amplitudeonly")
    return CalType::kDiagonalAmplitude;
  else if (modestr == "scalar" || modestr == "scalarcomplexgain" ||
           modestr == "scalarcomplex")
    return CalType::kScalar;
  else if (modestr == "scalaramplitude")
    return CalType::kScalarAmplitude;
  else if (modestr == "scalarphase")
    return CalType::kScalarPhase;
  else if (modestr == "tecandphase")
    return CalType::kTecAndPhase;
  else if (modestr == "tec")
    return CalType::kTec;
  else if (modestr == "tecscreen")
    return CalType::kTecScreen;
  else if (modestr == "fulljones")
    return CalType::kFullJones;
  else if (modestr == "rotation+diagonal")
    return CalType::kRotationAndDiagonal;
  else if (modestr == "rotation")
    return CalType::kRotation;
  throw std::runtime_error("Unknown mode: " + modestr);
}

std::string ToString(CalType caltype) {
  switch (caltype) {
    case CalType::kDiagonal:
      return "diagonal";
    case CalType::kScalar:
      return "scalarcomplexgain";
    case CalType::kFullJones:
      return "fulljones";
    case CalType::kDiagonalPhase:
      return "diagonalphase";
    case CalType::kScalarPhase:
      return "scalarphase";
    case CalType::kDiagonalAmplitude:
      return "diagonalamplitude";
    case CalType::kScalarAmplitude:
      return "scalaramplitude";
    case CalType::kTecAndPhase:
      return "tecandphase";
    case CalType::kTec:
      return "tec";
    case CalType::kTecScreen:
      return "tecscreen";
    case CalType::kRotation:
      return "rotation";
    case CalType::kRotationAndDiagonal:
      return "rotation+diagonal";
    default:
      throw std::runtime_error("Unknown caltype: " +
                               std::to_string(static_cast<int>(caltype)));
  }
}

}  // namespace base
}  // namespace dp3