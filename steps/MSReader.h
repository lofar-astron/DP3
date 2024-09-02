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

class MSReader : public InputStep {
 public:
  /// Construct the object for the given MS.
  /// Parameters are obtained from the parset using the given prefix.
  /// The missingData argument is for MultiMSReader.
  MSReader(const casacore::MeasurementSet& ms,
           const common::ParameterSet& parset, const std::string& prefix,
           bool missingData = false);

  /// Process the next data chunk.
  /// It returns false when at the end.
  bool process(std::unique_ptr<base::DPBuffer> buffer) override;

  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  /// Do nothing, since MSReader sets its info in the constructor.
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
  const casacore::Table& table() const override { return itsMS; }

  /// Get the slicer in the FLAG and DATA column.
  const casacore::Slicer& colSlicer() const { return itsColSlicer; }

  /// Get the rownrs for meta info of missing time slots.
  /// It uses the rows of the first time slot.
  const casacore::Vector<common::rownr_t>& getBaseRowNrs() const {
    return itsBaseRowNrs;
  }

  /// Get the name of the MS.
  std::string msName() const override;

  /// Get the baseline selection.
  const std::string& baselineSelection() const { return itsSelBL; }

  /// Is the data column missing?
  bool missingData() const { return itsMissingData; }

  /// Get the name of the main data column.
  const std::string& DataColumnName() const { return itsDataColName; }

  /// Get the names of the extra data columns that are also read.
  const std::vector<std::string>& ExtraDataColumnNames() const {
    return itsExtraDataColNames;
  }

  /// Get the weight column name.
  const std::string& WeightColumnName() const { return itsWeightColName; }

  /// Flags inf and NaN
  static void flagInfNaN(base::DPBuffer& buffer,
                         base::FlagCounter& flagCounter);

  /// Parse the expressions for the start channel and number of channels.
  /// The total number of channels can be used as 'nchan' in both expressions.
  static std::pair<unsigned int, unsigned int> ParseChannelSelection(
      const std::string& start_channel_string,
      const std::string& n_channels_string, unsigned int n_channels);

 protected:
  /// Default constructor.
  MSReader() = default;

 private:
  /// Prepare the access to the MS.
  /// Return the first and last time and the interval.
  void prepare(double& firstTime, double& lastTime, double& interval);

  /// Do the rest of the preparation.
  void prepare2(int spectralWindow, unsigned int start_channel,
                unsigned int n_channels);

  /// Skip the first times in the MS in case a start time was given.
  /// If needed, it sets itsFirstTime properly.
  void skipFirstTimes();

  /// Calculate the weights from the autocorrelations.
  void autoWeight(base::DPBuffer& buf);

  /// Read the weights at the given row numbers into the buffer.
  /// Note: the buffer must contain DATA if autoweighting is in effect.
  void getWeights(const casacore::RefRows& rowNrs, base::DPBuffer&);

 private:
  casacore::MeasurementSet itsMS;
  casacore::Table itsSelMS;  ///< possible selection of spw, baseline
  casacore::TableIterator itsIter;
  std::string itsDataColName{"DATA"};
  std::vector<std::string> itsExtraDataColNames{};
  std::string itsFlagColName{"FLAG"};
  std::string itsWeightColName{"WEIGHT_SPECTRUM"};
  std::string itsModelColName{"MODEL_DATA"};

 protected:
  std::string itsStartChanStr{"0"};  ///< startchan expression
  std::string itsNrChanStr{"0"};     ///< nchan expression
  std::string itsSelBL{};            ///< Baseline selection string
  bool itsNeedSort{false};           ///< sort needed on time,baseline?
  bool itsAutoWeight{false};         ///< calculate weights from autocorr?
  bool itsAutoWeightForce{false};    ///< always calculate weights?
  bool itsHasWeightSpectrum;
  bool itsUseFlags{true};
  bool itsUseAllChan{false};   ///< all channels (i.e. no slicer)?
  bool itsMissingData{false};  ///< allow missing data column?
  /// tolerance for time comparison
  ///
  /// Can be negative to insert flagged time slots before start.
  double itsTimeTolerance{1e-2};
  double itsTimeInterval{0.0};
  double itsFirstTime{0.0};
  double itsMaximumTime{0.0};  ///< the last/max. expected timeslot in the MS
  double itsNextTime{0.0};     ///< target time slot for the process() call
  double itsPrevMSTime{0.0};   ///< time of prev. iteration in search for target
  unsigned int itsNrRead{0};   ///< nr of time slots read from MS
  unsigned int itsNrInserted{0};  ///< nr of inserted time slots
  casacore::Slicer itsColSlicer;  ///< slice in corr,chan column
  casacore::Slicer itsArrSlicer;  ///< slice in corr,chan,bl array
  std::unique_ptr<base::UVWCalculator> itsUVWCalc;
  casacore::Vector<common::rownr_t>
      itsBaseRowNrs;  ///< rownrs for meta of missing times
  base::FlagCounter itsFlagCounter;
  common::NSTimer itsTimer;
};

}  // namespace steps
}  // namespace dp3

#endif
