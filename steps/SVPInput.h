// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Input step from streaming visibilities

#ifndef DP3_STEPS_SVPINPUT_H_
#define DP3_STEPS_SVPINPUT_H_

#include <steps/InputStep.h>

#include <dp3/steps/Step.h>
#include <dp3/base/DPBuffer.h>
#include <base/SubtableWriter.h>
#include <casacore/measures/Measures.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MCEpoch.h>
#include <common/ParameterSet.h>
#include <common/Timer.h>

#include <memory>
#include <string>
#include <vector>

namespace dp3 {
namespace steps {

class SVPInput : public InputStep {
  class MetaData {
   public:
    MetaData(){};
    size_t nr_channels_;
    int nr_polarizations_;
    size_t nr_antennas_;

    double start_time_;
    double end_time_;
    double integration_time_;
    size_t nr_times_;

    std::vector<casacore::MPosition> antenna_positions_;
    std::array<double, 9> axes_;
    std::vector<std::string> antenna_names_;
    std::vector<double> antenna_diameters_;
    double scale_factor_{1.0};

    double ra_, dec_;
    std::string source_name_;
    std::string frame_name_;
    std::string telescope_name_;

    std::vector<double> chan_freqs_;
    std::vector<double> chan_widths_;

   private:
  };

 public:
  SVPInput(const common::ParameterSet&, const std::string& prefix);

  void finish() override;
  bool process(std::unique_ptr<dp3::base::DPBuffer> buffer) override;
  // The following should be explicitly defined
  void updateInfo(const base::DPInfo&) override{};
  void show(std::ostream&) const override;

  std::string msName() const override;

 private:
  std::string name_;
  std::string temp_out_ms_;
  std::string receiver_socket_;
  const std::string kSocketName{"/tmp/svpsock0"};

  const std::string kDataColumnName{"DATA"};
  const std::string kWeightColumnName{"WEIGHT"};
  const std::string kFlagColumnName{"FLAG"};

  casacore::MDirection phase_direction_;
  std::unique_ptr<base::UVWCalculator> uvw_calculator_;

  MetaData metadata_;
  void InitializeInfo();
  void CreateInitialSubtables();
  void CreateAntennaTable(double time, base::SubtableWriter& writer);

  void CreateSpectralWindowTable(int nchannels,
                                 std::vector<double> const& chan_freqs,
                                 std::vector<double> const& chan_widths,
                                 base::SubtableWriter& writer);

  void CreateSourceTable(std::string& name, double start_time, double end_time,
                         base::SubtableWriter& writer);

  void CreateFieldTable(std::string& name, double time,
                        base::SubtableWriter& writer);

  void CreateObservationTable(double start_time, double end_time,
                              base::SubtableWriter& writer);

  void CreatePolarizationTable(base::SubtableWriter& writer,
                               const int nr_polarizations);

  int sock_fd_{-1};
  void ConnectServer(const char* path);
  void DisconnectServer() {
    if (sock_fd_ != -1) {
      close(sock_fd_);
    }
  }

  template <typename T>
  int ReceiveItem(T& message);
  int ReceiveBytes(char* message, const int message_size);
};

}  // namespace steps
}  // namespace dp3

#endif /* DP3_STEPS_SVPINPUT_H_ */
