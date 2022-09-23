// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_STEPS_OUTPUTSTEP_H_
#define DP3_STEPS_OUTPUTSTEP_H_

#include "Step.h"

namespace dp3 {
namespace steps {

/**
 * @brief Base class for output steps.
 */
class OutputStep : public Step {
 public:
  virtual ~OutputStep() {}

  /**
   * Set which fields the step should write.
   * @param fields A combination of fields. Non-writable fields are ignored.
   */
  void SetFieldsToWrite(const dp3::common::Fields& fields) {
    fields_to_write_ = fields;
  };

  /**
   * @return The fields the step should write.
   */
  const dp3::common::Fields& GetFieldsToWrite() const {
    return fields_to_write_;
  };

 private:
  /**
   * Determines which fields the step should write. Used by derived classes.
   */
  dp3::common::Fields fields_to_write_;
};

}  // namespace steps
}  // namespace dp3

#endif
