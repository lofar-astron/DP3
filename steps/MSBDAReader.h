// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Step for reading BDA data from an MS.
/// @author Lars Krombeen

#ifndef DPPP_MSBDAREADER_H
#define DPPP_MSBDAREADER_H

#include "InputStep.h"

#include "../base/BDABuffer.h"
#include "../base/UVWCalculator.h"
#include "../base/FlagCounter.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <casacore/tables/Tables/TableIter.h>
#include <casacore/tables/Tables/RefRows.h>
#include <casacore/casa/Arrays/Slicer.h>

#include <map>

namespace dp3 {
namespace common {
class ParameterSet;
}

namespace steps {
/// @brief Step for reading BDA data from an MS.

/// This class is a InputStep step reading the data from a MeasurementSet.
/// At the beginning it finds out the shape of the data; i.e., the
/// number of correlations, channels, baselines, and time slots.
/// It requires the data to be regularly shaped.
///
/// The object is constructed from the 'msin' keywords in the parset file.
/// Currently the following can be given:
/// <ul>
///  <li> msin: name of the MS
///  <li> msin.datacolumn: the data column to use [DATA]
///  <li> msin.weightcolumn: the weights column to use [WEIGHT_SPECTRUM or
///           WEIGHT]
/// </ul>
///
/// The data columns are handled in the following way:
/// <table>
///  <tr>
///   <td>TIME</td>
///   <td>The time slot center of the current data (in MJD seconds).
///       It is assumed that all data have the same interval, which is
///       used to find missing time slots.
///   </td>
///  </tr>
///  <tr>
///   <td>DATA</td>
///   <td>The visibility data as [ncorr,nchan,nbaseline]. Only the
///       part given by startchan and nchan is read. If a time slot is
///       inserted, all its data are zero.
///   </td>
///  </tr>
///  <tr>
///   <td>FLAG</td>
///   <td>The data flags as [ncorr,nchan,nbaseline] (True is bad).
///       They are read from the FLAG column. If a FLAG_ROW is set, all
///       flags for that baseline will be set. Also the flag of data
///       containing NaN or infinite numbers will be set.
///       All flags of an inserted time slot are set.
///   </td>
///  </tr>
///  <tr>
///   <td>WEIGHT</td>
///   <td>The data weights as [ncorr,nchan,nbaseline]. Column
///       WEIGHT_SPECTRUM is used if present and containing valid data,
///       otherwise column WEIGHT is used. The weights of an inserted
///       time slot are set to 0.
///       If autoweight is on, the autocorrelations are used to
///       calculate proper weights.
///   </td>
///  </tr>
///  <tr>
///   <td>UVW</td>
///   <td>The UVW coordinates in meters as [3,nbaseline].
///       They are calculated for a missing time slot.
///   </td>
///  </tr>
///  <tr>
///   <td>FULLRESFLAG</td>
///   <td>NOT YET IMPLEMENTED FOR BDA BUFFERS</td>
///   <td>For each baseline the LOFAR_FULL_RES_FLAG column is stored as
///       a uChar array with shape [orignchan/8, ntimeavg]. The bits
///       represent the flags. They are converted to a Bool array with shape
///       [orignchan, ntimeavg, nbaseline].
///       If column LOFAR_FULL_RES_FLAG is not present, the flags are used
///       and it is assumed that no averaging was done yet (thus ntimeavg=1
///       and orignchan=nchan).
///   </td>
///  </tr>
/// </table>

class MSBDAReader : public InputStep {
 public:
  /// Construct the object for the given MS.
  /// Parameters are obtained from the parset using the given prefix.
  MSBDAReader(const casacore::MeasurementSet& ms, const common::ParameterSet&,
              const std::string& prefix);

  virtual ~MSBDAReader();

  /// Process the next data chunk.
  /// It returns false when at the end.
  bool process(const base::DPBuffer&) override;

  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  /// Return which datatype this step outputs.
  MsType outputs() const override { return MsType::kBda; };

  /// Update the general info.
  void updateInfo(const base::DPInfo&) override;

  /// Add some data to the MeasurementSet written/updated.
  /// Do nothing.
  void addToMS(const std::string&) override{};

  /// Show the step parameters.
  void show(std::ostream&) const override;

  /// Show the timings.
  void showTimings(std::ostream&, double duration) const override;

  /// Get the name of the MS.
  std::string msName() const override;

  /// Tell if the visibility data are to be read. If set to true once,
  /// this will stay true.
  void setReadVisData(bool readVisData) override;

  /// Get the main MS table.
  const casacore::Table& table() const override { return ms_; }

  /// Get the selected spectral window.
  unsigned int spectralWindow() const override { return spw_; }

  /// Get the nr of averaged full resolution channels.
  unsigned int nchanAvgFullRes() const override { return 1; }

  /// Get the nr of averaged full resolution timeslots.
  unsigned int ntimeAvgFullRes() const override { return 1; }

  std::unique_ptr<everybeam::telescope::Telescope> GetTelescope(
      const everybeam::ElementResponseModel element_response_model,
      bool use_channel_frequency) const final override;

  /// Get the time information: cetroid of first and last time slots in the
  /// measurement set
  double firstTime() const override { return first_time_; }
  double lastTime() const override { return last_time_; }

 private:
  /// Reads the BDA subtables from an MS and stores the values that are required
  void FillInfoMetaData();

  const casacore::MeasurementSet ms_;
  std::string data_col_name_;
  std::string weight_col_name_;
  bool read_vis_data_;  ///< read visibility data?
  double last_ms_time_;
  double last_ms_interval_;
  double interval_;  ///< original interval of the MS
  bool is_interval_integer_;
  unsigned int spw_;    ///< spw (band) to use (<0 no select)
  unsigned int nread_;  ///< nr of time slots read from MS
  common::NSTimer timer_;
  std::size_t pool_size_;  ///< Pool size that will be used for the BDA buffers

  std::map<int, std::size_t>
      desc_id_to_nchan_;  ///< Maps DATA_DESC_ID to channel count.
  std::map<std::pair<int, int>, unsigned int>
      bl_to_id_;  ///< Maps a baseline(antenna pair) to a baseline index.

  double first_time_;
  double last_time_;
};

}  // namespace steps
}  // namespace dp3

#endif
