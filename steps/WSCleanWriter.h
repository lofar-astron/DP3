// WSCleanWriter.h: DP3 step writing a reordered MS
// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_STEPS_WSCLEAN_REORDER_WRITER_H_
#define DP3_STEPS_WSCLEAN_REORDER_WRITER_H_

#include "OutputStep.h"

#include <aocommon/logger.h>
#include <aocommon/io/serialstreamfwd.h>
#include <aocommon/polarization.h>

#include <schaapcommon/reordering/reorderedfilewriter.h>

#include "../base/FlagCounter.h"
#include "../common/ParameterSet.h"
#include "../common/Timer.h"

namespace dp3 {

namespace common {
class ParameterSet;
}

namespace steps {

class WSCleanWriter : public OutputStep {
 public:
  const common::Fields kRequiredFields =
      kDataField | kFlagsField | kWeightsField | kUvwField;

  explicit WSCleanWriter(const common::ParameterSet& parset,
                         const std::string& prefix);
  ~WSCleanWriter();

  common::Fields getRequiredFields() const override { return kRequiredFields; }

  bool process(std::unique_ptr<base::DPBuffer> buffer) override;

  void finish() override;

  void show(std::ostream& os) const override;

  void updateInfo(const base::DPInfo& info_obj) override;

  void showTimings(std::ostream& os, double duration) const override;

 private:
  void StartReorder();
  void ReorderBuffer(dp3::base::DPBuffer& buffer);
  void FinishReorder();

  std::string name_;
  common::ParameterSet parset_;
  std::string out_name_;

  std::string temporary_directory_;
  std::set<aocommon::PolarizationEnum> pols_out_;

  uint32_t data_desc_id_;

  dp3::base::DPInfo dp_info_;

  std::unique_ptr<schaapcommon::reordering::ReorderedFileWriter> writer_;

  common::NSTimer timer_;
  common::NSTimer writer_timer_;
};
}  // namespace steps
}  // namespace dp3

#endif
