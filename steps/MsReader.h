// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DP3 step reading from an MS
/// @author Ger van Diepen

#ifndef DP3_STEPS_MSREADER_H_
#define DP3_STEPS_MSREADER_H_

#include <casacore/casa/Arrays/Slicer.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/tables/Tables/TableIter.h>

#include <dp3/base/DPBuffer.h>

#include "InputStep.h"

#include "../base/UVWCalculator.h"
#include "../base/FlagCounter.h"

namespace dp3 {
namespace steps {
/// @brief DP3 step reading from an MS

/// This class is an InputStep step reading the data from a MeasurementSet.
/// At the beginning it finds out the shape of the data; i.e., the
/// number of correlations, channels, baselines, and time slots.
/// It requires the data to be regularly shaped.
///
/// The object is constructed from the 'msin' keywords in the parset file.
/// Currently the following can be given:
/// <ul>
///  <li> msin: name of the MS
///  <li> msin.autoweight: calculate weights from autocorrelations? [no]
///  <li> msin.startchan: first channel to use [0]
///  <li> msin.nchan: number of channels to use [all]
///  <li> msin.useflag: use the existing flags? [yes]
///  <li> msin.datacolumn: the data column to use (measured visibilities) [DATA]
///  <li> msin.extradatacolumns: the extra data columns to use. Could be used to
///           load model visibilities when stored in the MS as well [empty]
///  <li> msin.weightcolumn: the weights column to use [WEIGHT_SPECTRUM or
///           WEIGHT]
///  <li> msin.starttime: first time to use [first time in MS]
///  <li> msin.endtime: last time to use [last time in MS]
/// </ul>
///
/// If a time slot is missing, it inserts flagged data set to zero, i.e. it
/// inserts dummy data that is marked invalid.
/// Missing time slots can also be detected at the beginning or end of the
/// MS by giving the correct starttime and endtime.
/// The correct UVW coordinates are calculated for inserted time slots.
///
/// The process function reads the fields needed by the chain of steps which can
/// be obtained calling getFieldsToRead(). In this way only the necessary data
/// is kept in memory.
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
/// </table>

class MsReader : public InputStep {
 public:
  /// Construct the object for the given MS.
  /// Parameters are obtained from the parset using the given prefix.
  /// The allow_missing_data argument is for MultiMsReader.
  MsReader(const casacore::MeasurementSet& ms,
           const common::ParameterSet& parset, const std::string& prefix,
           bool allow_missing_data = false);

  /// Process the next data chunk.
  /// It returns false when at the end.
  bool process(std::unique_ptr<base::DPBuffer> buffer) override;

  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  /// Do nothing, since MsReader sets its info in the constructor.
  void updateInfo(const base::DPInfo&) override {}

  /// Add some data to the MeasurementSet written/updated.
  /// Do nothing.
  void addToMS(const std::string&) override{};

  /// Show the step parameters.
  void show(std::ostream&) const override;

  /// If needed, show the flag counts.
  void showCounts(std::ostream&) const override;

  /// Show the timings.
  void showTimings(std::ostream&, double duration) const override;

  /// Read the UVW at the given row numbers into the buffer.
  void getUVW(const casacore::RefRows& rowNrs, double time, base::DPBuffer&);

  /// Get the main MS table.
  const casacore::Table& table() const override { return ms_; }

  /// Get the name of the MS.
  std::string msName() const override;

  /// Get the baseline selection.
  const std::string& BaselineSelection() const { return baseline_selection_; }

  /// Is the data column missing?
  bool MissingData() const { return missing_data_; }

  /// Get the names of the extra data columns that are also read.
  const std::vector<std::string>& ExtraDataColumnNames() const {
    return extra_data_column_names_;
  }

  bool AutoWeight() const { return auto_weight_; }

  /// Flags infinite and Not-a-Number values.
  static void FlagInfinityNan(base::DPBuffer& buffer,
                              base::FlagCounter& flagCounter);

  /// Parse the expressions for the start channel and number of channels.
  /// The total number of channels can be used as 'nchan' in both expressions.
  static std::pair<unsigned int, unsigned int> ParseChannelSelection(
      const std::string& start_channel_string,
      const std::string& n_channels_string, unsigned int n_channels);

 private:
  /// Applies baseline_selection_ to selection_ms_, if needed.
  void SelectBaselines();

  /// Performs MS column related initialization.
  /// - Sets has_weight_spectrum_ and missing_data_.
  /// - Stores MS and column names in infoOut().
  void InitializeColumns(bool allow_missing_data,
                         const std::string& data_column_name,
                         const std::string& flag_column_name,
                         const std::string& weight_column_name);

  /// Initializes ms_iterator_ for selection_ms_.
  /// If needed, creates the iterator an a sorted MS.
  /// Reads the first time, last time and interval from the MS and stores them
  /// in infoOut().
  /// Also handles auto-weighting settings.
  /// Needs to happen after SelectBaselines(), which may update selection_ms_.
  /// @param force_auto_weight If true, always enable auto weighting.
  void InitializeIterator(bool force_auto_weight);

  /// Reads antenna information from the MS and stores it in infoOut().
  /// @param table The table for the MS with selection and sorting applied,
  ///        as initialized by InitializeIterator().
  void ReadAntennas(const casacore::Table& table);

  /// Reads array information from the MS and stores it in infoOut().
  /// Creates the UVW calculator using the array information.
  /// Needs to happen after ReadAntennas(), since it uses the antenna positions.
  void ReadArrayInformation();

  /// Processes time selection settings.
  /// Updates the first time and last time in infoOut() if needed.
  /// Needs to happen after InitializeIterator(), which initially stores the MS
  /// times into infoOut().
  void ParseTimeSelection(const common::ParameterSet& parset,
                          const std::string& prefix);

  /// Skip the first times in the MS in case a start time was given.
  /// Update first_time if needed.
  void SkipFirstTimes(double& first_time, double interval);

  /// Do all channel-related initialization.
  void InitializeChannels(int spectralWindow);

  /// Read frequency, width, resolution, effective bandwidth and
  /// reference frequency from the MS and call infoOut().setChannels().
  void ReadChannelProperties(int spectralWindow);

  /// Read polarizations from the MS and store in infoOut().
  void ReadPolarizations(int spectralWindow);

  /// Calculate the weights from the autocorrelations.
  void AutoWeight(base::DPBuffer& buf);

  /// Read the weights at the given row numbers into the buffer.
  /// Note: the buffer must contain DATA if autoweighting is in effect.
  void GetWeights(const casacore::RefRows& rowNrs, base::DPBuffer&);

  casacore::MeasurementSet ms_;
  casacore::Table selection_ms_;  ///< possible selection of spw, baseline
  casacore::TableIterator ms_iterator_;
  std::vector<std::string> extra_data_column_names_;
  std::string start_channel_expression_{"0"};
  std::string n_channels_expression_{"0"};
  std::string baseline_selection_{};
  bool sort_{false};         ///< sort needed on time,baseline?
  bool auto_weight_{false};  ///< calculate weights from autocorr?
  bool has_weight_spectrum_{false};
  bool use_flags_{true};
  bool use_all_channels_{false};  ///< all channels (i.e. no slicer)?
  bool missing_data_{false};      ///< is the data column missing?
  /// tolerance for time comparison
  ///
  /// Can be negative to insert flagged time slots before start.
  double time_tolerance_{1.0e-2};
  /// Correction on first time in case times are irregularly spaced.
  double time_correction_{0.0};
  double next_time_{0.0};      ///< target time slot for the process() call
  double previous_time_{0.0};  ///< time of prev. iteration in search for target
  unsigned int n_read_{0};     ///< nr of time slots read from MS
  unsigned int n_inserted_{0};      ///< nr of inserted time slots
  casacore::Slicer column_slicer_;  ///< slice in corr,chan column
  casacore::Slicer array_slicer_;   ///< slice in corr,chan,bl array
  std::unique_ptr<base::UVWCalculator> uvw_calculator_;
  base::FlagCounter flag_counter_;
  common::NSTimer timer_;
};

}  // namespace steps
}  // namespace dp3

#endif
