// ParameterRecord.h: A record of parameter values
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef LOFAR_COMMON_PARAMETERRECORD_H
#define LOFAR_COMMON_PARAMETERRECORD_H

#include "ParameterSet.h"

namespace dp3 {
namespace common {

/// \brief A record of parameter values
class ParameterRecord : public ParameterSet {
 public:
  /// Define the iterators for this class.
  typedef ParameterSet::iterator iterator;
  typedef ParameterSet::const_iterator const_iterator;

  /// Default constructor creates empty record.
  ParameterRecord() {}

  /// Try to get a value from the record or from a nested record.
  bool getRecursive(const std::string& key, ParameterValue& value) const;

  /// Put to ostream.
  friend std::ostream& operator<<(std::ostream& os, const ParameterRecord&);
};

}  // namespace common
}  // namespace dp3

#endif
