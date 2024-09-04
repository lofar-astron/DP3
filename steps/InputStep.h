// InputStep.h: Abstract base class for a Step generating input
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Abstract base class for a Step generating input
/// @author Ger van Diepen

#ifndef DP3_STEPS_INPUTSTEP_H_
#define DP3_STEPS_INPUTSTEP_H_

#include <memory>

#include <dp3/steps/Step.h>

#include "../base/UVWCalculator.h"
#include "../base/FlagCounter.h"
#include "../common/ParameterSet.h"
#include "../common/Timer.h"

namespace casacore {
class MeasurementSet;
class RefRows;
}  // namespace casacore

namespace dp3 {
namespace steps {

/// @brief Abstract base class for a Step generating input

/// This class is the abstract base class for a Step object that
/// handles the input. A concrete example is MsReader that reads the
/// data from a MeasurementSet. However, it is also possible to have
/// input steps generating data on the fly as done in test programs
/// like tAverager.cc.

class InputStep : public Step {
 public:
  ~InputStep() override;

  common::Fields getRequiredFields() const override { return {}; }

  common::Fields getProvidedFields() const override { return fields_to_read_; }

  /// Get the MS name.
  /// The default implementation returns an empty string.
  virtual std::string msName() const;

  /// Set which fields must be read.
  virtual void setFieldsToRead(const dp3::common::Fields& fields) {
    fields_to_read_ = fields;
  };

  /// Get which fields must be read.
  const dp3::common::Fields& getFieldsToRead() const {
    return fields_to_read_;
  };

  /// Get the main MS table.
  virtual const casacore::Table& table() const;

  /// Check if a measurement set contains Baseline Dependent Averaged data.
  /// @param ms A casacore measurement set.
  /// @return true if the measurement set has BDA data, false if it is regular.
  static bool HasBda(const casacore::MeasurementSet& ms);

  /// Creates a (multi) MS reader.
  /// If it receives a single input MS name, it will create either a regular
  /// MsReader step or a MSBDAReader step depending on the contents of the MS.
  /// If it receives multiple input MS names, it will create a MultiMsReader
  /// step. In this case, BDA data is not supported (yet).
  static std::unique_ptr<InputStep> CreateReader(const common::ParameterSet&);

 private:
  /// This variable is used by the inputStep's derived classes to determine
  /// which fields must be read.
  dp3::common::Fields fields_to_read_;
};

}  // namespace steps
}  // namespace dp3

#endif
