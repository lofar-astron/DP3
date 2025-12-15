// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_DDECAL_CONSTRAINT_RESULT_H_
#define DP3_DDECAL_CONSTRAINT_RESULT_H_

#include <string>
#include <vector>

namespace dp3::ddecal {

/**
 * Holds the information that is produced by a constraint and that
 * should be written to the h5 solution file. Examples are
 * FaradayConstraint and TECConstraint, which store their found faraday
 * or tec values in a ConstraintResult object.
 */
struct ConstraintResult {
  /// Both vals and weights have the dimensions described in dims and axes.
  std::vector<double> vals;
  std::vector<double> weights;
  /// Comma-separated string with axis names, fastest varying last.
  std::string axes;
  std::vector<size_t> dims;
  std::string name;
};

}  // namespace dp3::ddecal

#endif
