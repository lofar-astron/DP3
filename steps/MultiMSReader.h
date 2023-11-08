// MultiMSReader.h: DP3 step reading from multiple MSs
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DP3 step reading from multiple MSs
/// @author Ger van Diepen

#ifndef DP3_STEPS_MULTIMSREADER_H_
#define DP3_STEPS_MULTIMSREADER_H_

#include <casacore/tables/Tables/TableIter.h>
#include <casacore/tables/Tables/RefRows.h>
#include <casacore/casa/Arrays/Slicer.h>

#include "../base/UVWCalculator.h"
#include "../base/FlagCounter.h"
#include "../common/ParameterSet.h"

#include "MSReader.h"
#include "ResultStep.h"

namespace dp3 {
namespace steps {
/// @brief DP3 step reading from multiple MSs

/// This class is a InputStep step reading the data from a MeasurementSet.
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
///  <li> msin.datacolumn: the data column to use [DATA]
///  <li> msin.starttime: first time to use [first time in MS]
///  <li> msin.endtime: last time to use [last time in MS]
/// </ul>
///
/// If a time slot is missing, it is inserted with flagged data set to zero.
/// Missing time slots can also be detected at the beginning or end of the
/// MS by giving the correct starttime and endtime.
/// The correct UVW coordinates are calculated for inserted time slots.
///
/// The process function reads the fields needed by the steps' chain which can
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

class MultiMSReader final : public MSReader {
 public:
  /// Construct the object for the given MS.
  /// Parameters are obtained from the parset using the given prefix.
  MultiMSReader(const std::vector<string>& msNames, const common::ParameterSet&,
                const string& prefix);

  ~MultiMSReader() override;

  /// Process the next data chunk.
  /// It returns false when at the end.
  bool process(std::unique_ptr<base::DPBuffer> buffer) override;

  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  /// Update the general info (by initializing it).
  void updateInfo(const base::DPInfo&) override;

  /// Show the step parameters.
  void show(std::ostream&) const override;

  /// If needed, show the flag counts.
  void showCounts(std::ostream&) const override;

  /// Show the timings.
  void showTimings(std::ostream&, double duration) const override;

  /// Set which fields must be read.
  void setFieldsToRead(const dp3::common::Fields& fields) override;

  /// Get the name of the first MS.
  std::string msName() const override;

 private:
  /// Handle the info for all bands.
  void handleBands();

  /// Sort the bands (MSs) inorder of frequency.
  void sortBands();

  /// Fill the band info where some MSs are missing.
  void fillBands();

  /// Reads the weights into 'buffer'
  void getWeights(std::unique_ptr<base::DPBuffer>& buffer);

  bool itsOrderMS;  ///< sort multi MS in order of freq?
  int itsFirst;     ///< first valid MSReader (<0 = none)
  int itsNMissing;  ///< nr of missing MSs
  struct Reader {
    std::string name;
    std::shared_ptr<MSReader> ms_reader;
    std::shared_ptr<ResultStep> result;
  };
  std::vector<Reader> readers_;
  unsigned int itsFillNChan;  ///< nr of chans for missing MSs
  base::FlagCounter itsFlagCounter;
};

}  // namespace steps
}  // namespace dp3

#endif
