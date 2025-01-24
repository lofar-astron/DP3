// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// This class is responsible to create ancillary tables such as
/// ANTENNA or SPECTRAL_WINDOW
/// In the future this could be part of the DP3 writer step
/// it was a choice to keep it separate to avoid changing core
/// DP3 functionalities

/// Most of this code was taken from https://git.astron.nl/RD/aartfaac-tools.git
/// and adapted to fit into the DP3 codebase.

#ifndef DP3_BASE_SUBTABLEWRITER_H_
#define DP3_BASE_SUBTABLEWRITER_H_

#include <array>
#include <complex>
#include <string>
#include <vector>

#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/Containers/Record.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/DataMan/IncrementalStMan.h>
#include <casacore/tables/Tables/ArrColDesc.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScaColDesc.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/SetupNewTab.h>
#include <casacore/tables/Tables/TableRecord.h>

namespace dp3::base {

class SubtableWriter {
 public:
  struct AntennaInfo {
    std::string name, station;
    std::string type, mount;
    double x, y, z;
    double diameter;
    bool flag;
  };

  struct ChannelInfo {
    double channel_frequency;
    double channel_width;
    double effective_bandwidth;
    double resolution;
  };

  struct SourceInfo {
    int source_id;
    double time;
    double interval;
    int spectral_window_id;
    int num_lines;
    std::string name;
    int calibration_group;
    std::string code;
    double ra;
    double dec;
    double proper_motion[2];
  };

  struct FieldInfo {
    std::string name;
    std::string code;
    double time;
    int num_poly;
    double delay_direction_ra;
    double delay_direction_dec;
    double phase_direction_ra;
    double phase_direction_dec;
    double reference_direction_ra;
    double reference_direction_dec;
    int source_id;
    bool flag_row;
  };

  struct ObservationInfo {
    std::string telescope_name;
    double start_time;
    double end_time;
    std::string observer;
    std::string schedule_type;
    std::string project;
    double release_date;
    bool flag_row;
    std::string antenna_type;
    // RCU -> ReCeiver Unit
    int rcu_mode;
    int flag_window_size;
  };

  SubtableWriter(){};
  SubtableWriter(std::string path, const int nr_channels = 1);

  ~SubtableWriter(){};
  void WriteBandInfo(const std::string &name,
                     const std::vector<ChannelInfo> &channels,
                     double reference_frequency, double total_bandwidth,
                     bool flag_row);
  void WriteAntennas(const std::vector<AntennaInfo> &antennas,
                     const std::array<double, 9> &coordinate_axes, double time);
  void WriteLinearPolarizations(bool flagRow, const int n_pol = 4);
  void WriteSource(const SourceInfo &source);
  void WriteField(const FieldInfo &field);
  void WriteObservation(const ObservationInfo &observation);
  void WriteHistoryItem(const std::string &commandLine,
                        const std::string &application,
                        const std::vector<std::string> &params);
  const std::string &GetPath() { return path_; }

 private:
  void WriteDataDescEntry(size_t spectral_window_id, size_t polarization_id,
                          bool flag_row);
  void WriteFeedEntries(const std::vector<AntennaInfo> &antennas, double time);
  std::string path_;
  casacore::MeasurementSet ms_;
};
}  // namespace dp3::base
#endif  // DP3_BASE_SUBTABLEWRITER_H_
