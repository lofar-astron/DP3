// WSCleanWriter.h: DP3 step writing a reordered MS
// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_STEPS_WSCLEAN_REORDER_WRITER_H_
#define DP3_STEPS_WSCLEAN_REORDER_WRITER_H_

#include "OutputStep.h"

#include <cstddef>
#include <fstream>

#include <aocommon/logger.h>
#include <aocommon/io/serialstreamfwd.h>
#include <aocommon/polarization.h>

#include "../base/FlagCounter.h"
#include "../common/ParameterSet.h"
#include "../common/Timer.h"

namespace dp3 {

namespace common {
class ParameterSet;
}

namespace steps {

struct ReorderFile {
  std::unique_ptr<std::ofstream> data;
  std::unique_ptr<std::ofstream> weight;
};

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
  aocommon::PolarizationEnum polarization_;
  size_t nr_polarizations_;

  double start_time_;
  uint64_t channel_count_;
  uint64_t channel_start_;
  uint32_t data_desc_id_;

  size_t selected_row_count_;

  dp3::base::DPInfo dp_info_;

  std::unique_ptr<std::ofstream> meta_file_ptr_;
  std::vector<ReorderFile> files_;

  common::NSTimer timer_;
  common::NSTimer writer_timer_;
};
}  // namespace steps
}  // namespace dp3

#endif
