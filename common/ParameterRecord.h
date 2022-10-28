// ParameterRecord.h: A record of parameter values
//
// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_COMMON_PARAMETERRECORD_H_
#define DP3_COMMON_PARAMETERRECORD_H_

#include "ParameterSet.h"

namespace dp3 {
namespace common {

/// \brief A record of parameter values.
/// The only difference with a ParameterSet is the output operator.
class ParameterRecord : public ParameterSet {
 public:
  /// Put to ostream.
  friend std::ostream& operator<<(std::ostream& os, const ParameterRecord&);
};

}  // namespace common
}  // namespace dp3

#endif
