// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_STEPS_AARTFAAC_READER_H_
#define DP3_STEPS_AARTFAAC_READER_H_

#include "aartfaac/PacketHeader.h"
#include "base/DPBuffer.h"
#include "base/RcuMode.h"
#include "base/SubtableWriter.h"
#include "steps/InputStep.h"
#include "steps/Step.h"

#include <casacore/measures/Measures.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MCEpoch.h>
#include <common/ParameterSet.h>
#include <common/Timer.h>

#include <memory>
#include <string>
#include <vector>

namespace dp3::aartfaac {
class File;
}  // namespace dp3::aartfaac

namespace dp3::steps {

/// @brief DP3 step class that inputs data from the AARTFAAC correlator.
/// This class is meant to stream directly the data from the AARTFAAC
/// correlator in DP3
class AartfaacReader : public InputStep {
 public:
  /// Construct the object.
  /// Parameters are obtained from the parset using the given prefix.
  AartfaacReader(const common::ParameterSet&, const std::string& prefix);
  ~AartfaacReader();

  /// Finish the processing of this step and subsequent steps.
  void finish() override;
  bool process(std::unique_ptr<dp3::base::DPBuffer> buffer) override;
  void updateInfo(const base::DPInfo&) override;
  /// Show the step parameters.
  void show(std::ostream&) const override;

  std::string msName() const override;

 private:
  base::RcuMode mode_ = base::RcuMode::Unused;
  std::string antenna_configuration_;
  std::string out_name_;
  std::string stream_type_;
  std::string filename_;
  std::string data_col_name_ = "DATA";
  std::string flag_col_name_ = "FLAG";
  std::string weight_col_name_ = "WEIGHT";
  common::NSTimer timer_;
  std::unique_ptr<dp3::aartfaac::File> input_;

  std::vector<casacore::MPosition> antenna_positions_;
  std::array<double, 9> axes_;
  casacore::MDirection phase_direction_;
  std::vector<std::string> antenna_names_;
  std::vector<double> antenna_diameters_;
  dp3::aartfaac::PacketHeader metadata_;
  std::unique_ptr<base::UVWCalculator> uvw_calculator_;
  double start_time_ = 0.0;
  double end_time_ = 0.0;
  double integration_time_ = 0.0;
  int current_time_ = 0;
  int n_times_ = 0;

  void CreateInitialSubtables();
  void InitializeInfo();
  void OpenFile();

  void ApplyDelayCorrection(dp3::base::DPBuffer& buffer, size_t baseline);

  void ComputeWeights(dp3::base::DPBuffer& buffer, size_t baseline,
                      float exposure);

  void CreateAntennaTable(double time, base::SubtableWriter& writer);

  void CreateSpectralWindowTable(int nchannels, double start_frequency,
                                 double channel_width,
                                 base::SubtableWriter& writer);

  void CreateSourceTable(double start_time, double end_time,
                         base::SubtableWriter& writer);

  void CreateFieldTable(double time, base::SubtableWriter& writer);

  void CreateObservationTable(double start_time, double end_time,
                              base::SubtableWriter& writer);

  void CreatePolarizationTable(base::SubtableWriter& writer);
};

}  // namespace dp3::steps

#endif  // DP3_STEPS_AARTFAAC_READER_H_
