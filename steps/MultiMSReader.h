// MultiMSReader.h: DPPP step reading from multiple MSs
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DPPP step reading from multiple MSs
/// @author Ger van Diepen

#ifndef DPPP_MULTIMSREADER_H
#define DPPP_MULTIMSREADER_H

#include "MSReader.h"

#include "../base/DPBuffer.h"
#include "../base/UVWCalculator.h"
#include "../base/FlagCounter.h"

#include <casacore/tables/Tables/TableIter.h>
#include <casacore/tables/Tables/RefRows.h>
#include <casacore/casa/Arrays/Slicer.h>

namespace dp3 {
namespace common {
class ParameterSet;
}

namespace steps {
/// @brief DPPP step reading from multiple MSs

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
/// The process function only reads the data and flags to avoid that
/// too much data is kept in memory.
/// Other columns (like WEIGHT, UVW) can be read when needed by using the
/// appropriate InputStep::fetch function.
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

class MultiMSReader final : public MSReader {
 public:
  /// Construct the object for the given MS.
  /// Parameters are obtained from the parset using the given prefix.
  MultiMSReader(const std::vector<string>& msNames, const common::ParameterSet&,
                const string& prefix);

  ~MultiMSReader() override;

  /// Process the next data chunk.
  /// It returns false when at the end.
  bool process(const base::DPBuffer&) override;

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

  /// Read the UVW at the given row numbers.
  void getUVW(const casacore::RefRows& rowNrs, double time,
              base::DPBuffer& buf) override;

  /// Read the weights at the given row numbers.
  void getWeights(const casacore::RefRows& rowNrs,
                  base::DPBuffer& buf) override;

  /// Read the FullRes flags (LOFAR_FULL_RES_FLAG) at the given row numbers.
  /// It returns a 3-dim array [norigchan, ntimeavg, nbaseline].
  /// If undefined, false is returned.
  bool getFullResFlags(const casacore::RefRows& rowNrs,
                       base::DPBuffer& buf) override;

  /// Tell if the visibility data are to be read.
  void setReadVisData(bool readVisData) override;

  /// Get the name of the first MS.
  std::string msName() const override;

 private:
  /// Handle the info for all bands.
  void handleBands();

  /// Sort the bands (MSs) inorder of frequency.
  void sortBands();

  /// Fill the band info where some MSs are missing.
  void fillBands();

  bool itsOrderMS;  ///< sort multi MS in order of freq?
  int itsFirst;     ///< first valid MSReader (<0 = none)
  int itsNMissing;  ///< nr of missing MSs
  std::vector<string> itsMSNames;
  std::vector<std::shared_ptr<MSReader>> itsReaders;
  std::vector<base::DPBuffer> itsBuffers;
  unsigned int itsFillNChan;  ///< nr of chans for missing MSs
  base::FlagCounter itsFlagCounter;
  bool itsRegularChannels;  /// Are resulting channels regularly spaced
};

}  // namespace steps
}  // namespace dp3

#endif
