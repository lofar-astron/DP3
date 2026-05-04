// Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "CalType.h"

#include <stdexcept>

namespace dp3::base {

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
  else if (modestr == "tec")
    return CalType::kTec;
  else if (modestr == "tec+delay")
    return CalType::kTecAndDelay;
  // tecandphase is the old name, and is still supported for compatibility
  else if (modestr == "tec+phase" || modestr == "tecandphase")
    return CalType::kTecAndPhase;
  else if (modestr == "tec+phase+delay")
    return CalType::kTecPhaseAndDelay;
  else if (modestr == "tecscreen")
    return CalType::kTecScreen;
  else if (modestr == "fulljones")
    return CalType::kFullJones;
  else if (modestr == "rotation+diagonal")
    return CalType::kRotationAndDiagonal;
  else if (modestr == "rotation")
    return CalType::kRotation;
  else if (modestr == "faradayrotation")
    return CalType::kFaradayRotation;
  else if (modestr == "leakage")
    return CalType::kLeakage;
  else if (modestr == "leakageamplitude")
    return CalType::kLeakageAmplitude;
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
    case CalType::kTec:
      return "tec";
    case CalType::kTecAndDelay:
      return "tec+delay";
    case CalType::kTecAndPhase:
      return "tec+phase";
    case CalType::kTecPhaseAndDelay:
      return "tec+phase+delay";
    case CalType::kTecScreen:
      return "tecscreen";
    case CalType::kRotation:
      return "rotation";
    case CalType::kRotationAndDiagonal:
      return "rotation+diagonal";
    case CalType::kFaradayRotation:
      return "faraday";
    case CalType::kLeakage:
      return "leakage";
    case CalType::kLeakageAmplitude:
      return "leakageamplitude";
    default:
      throw std::runtime_error("Unknown caltype: " +
                               std::to_string(static_cast<int>(caltype)));
  }
}

}  // namespace dp3::base
