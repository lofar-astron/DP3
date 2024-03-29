// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_STEPS_OUTPUTSTEP_H_
#define DP3_STEPS_OUTPUTSTEP_H_

#include <dp3/steps/Step.h>

namespace dp3 {
namespace steps {

/**
 * @brief Base class for output steps.
 */
class OutputStep : public Step {
 public:
  /**
   * Common override for all output steps. Since they should output all data,
   * they never provide any new fields.
   */
  common::Fields getProvidedFields() const override { return {}; }

  /**
   * Set which fields the step should write.
   * @param fields A combination of fields. Non-writable fields are ignored.
   */
  virtual void SetFieldsToWrite(const dp3::common::Fields& fields) {
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
