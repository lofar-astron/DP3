// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// This code was taken from
/// https://git.astron.nl/RD/aartfaac-tools/-/raw/master/common/aartfaac/aartfaacmode.h?ref_type=heads
/// and adapted to fit into DP3 codebase.
///
/// The Receiver Unit Mode, hereafter RCU mode, is the mode in which the
/// receiver units are set in a LOFAR station. The mode determines the amount
/// of frequency ranges the RCU is sensitive towards. The RcuMode class
/// provides the the frequency, bandwidth, antenna type, and a descriptive
/// string of the RCU mode given a mode code.
/// A reference for this modes can be found in
/// https://science.astron.nl/telescopes/lofar/lofar-system-overview/technical-specification/frequency-subband-selection-and-rfi/.

#ifndef DP3_BASE_RCU_MODE_H
#define DP3_BASE_RCU_MODE_H

#include <string>
#include <stdexcept>

namespace dp3::base {

class RcuMode {
 public:
  enum Mode {
    Unused = 0,         // 0 = Unused
    LBAOuter10_90 = 1,  // 1 = LBA_OUTER, 10-90 MHz Analog filter
    LBAOuter30_90 = 2,  // Mode 2 = LBA_OUTER, 30-90 MHz Analog filter
    LBAInner10_90 = 3,  // Mode 3 = LBA_INNER, 10-90 MHz Analog filter
    LBAInner30_90 = 4,  // Mode 4 = LBA_INNER, 30-90 MHz Analog filter
    HBA110_190 = 5,     // Mode 5 = HBA, 110-190MHz Analog filter
    HBA170_230 = 6,     // Mode 6 = HBA, 170-230MHz Analog filter
    HBA210_270 = 7      // Mode 7 = HBA, 210-270MHz Analog filter
  } mode;

  RcuMode(const Mode& m) : mode(m){};
  RcuMode(){};

  std::string ToString() const {
    switch (mode) {
      default:
      case Unused:
        return "unused";
      case LBAOuter10_90:
        return "LBA_OUTER 10-90 MHz";
      case LBAOuter30_90:
        return "LBA_OUTER 30-90 MHz";
      case LBAInner10_90:
        return "LBA_INNER 10-90 MHz";
      case LBAInner30_90:
        return "LBA_INNER 30-90 MHz";
      case HBA110_190:
        return "HBA 110-190 MHz";
      case HBA170_230:
        return "HBA 170-230 MHz";
      case HBA210_270:
        return "HBA 210-270 MHz";
    }
  }

  std::string AntennaType() {
    std::string antenna_type;
    switch (mode) {
      case RcuMode::LBAInner10_90:
      case RcuMode::LBAInner30_90:
      case RcuMode::LBAOuter10_90:
      case RcuMode::LBAOuter30_90:
        antenna_type = "LBA";
        break;
      case RcuMode::HBA110_190:
      case RcuMode::HBA170_230:
      case RcuMode::HBA210_270:
        antenna_type = "HBA";
        break;
      default:
        antenna_type = "?";
        break;
    }
    return antenna_type;
  }
  bool operator==(const Mode& _mode) const { return _mode == mode; }

  static RcuMode FromNumber(const int& modeNumber) {
    return {static_cast<Mode>(modeNumber)};
  }

  double Bandwidth() const {
    switch (mode) {
      case RcuMode::LBAInner10_90:  // 200 MHz clock, Nyquist zone 1
      case RcuMode::LBAInner30_90:
      case RcuMode::LBAOuter10_90:
      case RcuMode::LBAOuter30_90:
      case RcuMode::HBA110_190:  // 200 MHz clock, Nyquist zone 2
      case RcuMode::HBA210_270:  // 200 MHz clock, Nyquist zone 3
        return 195312.5;         // 1/1024 x nu_{clock}
      case RcuMode::HBA170_230:  // 160 MHz clock, Nyquist zone 3
        return 156250.0;
      default:
        throw std::runtime_error(
            "Don't know how to handle this mode: not implemented yet");
    }
  }

  double CentralFrequency() const {
    switch (mode) {
      case RcuMode::LBAOuter30_90:
      case RcuMode::LBAInner30_90:
        return 60.0;
      case RcuMode::LBAOuter10_90:
      case RcuMode::LBAInner10_90:
        return 50.0;
      case RcuMode::HBA110_190:
        return 150.0;
      case RcuMode::HBA210_270:
        return 240.0;
      case RcuMode::HBA170_230:
        return 200.0;
      default:
        throw std::runtime_error(
            "Don't know how to handle this mode: not implemented yet");
    }
  }

  double FrequencyOffset() const {
    switch (mode) {
      case RcuMode::LBAInner10_90:  // 200 MHz clock, Nyquist zone 1
      case RcuMode::LBAInner30_90:
      case RcuMode::LBAOuter10_90:
      case RcuMode::LBAOuter30_90:
        return 0.0;
      case RcuMode::HBA110_190:  // 200 MHz clock, Nyquist zone 2
        return 100.0e6;
      case RcuMode::HBA170_230:  // 160 MHz clock, Nyquist zone 3
        return 160.0e6;
      case RcuMode::HBA210_270:  // 200 MHz clock, Nyquist zone 3
        return 200.0e6;
      default:
        throw std::runtime_error(
            "Don't know how to handle this mode: not implemented yet");
    }
  }
};
}  // namespace dp3::base
#endif  // DP3_BASE_RCU_MODE_H
