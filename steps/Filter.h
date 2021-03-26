// Filter.h: DPPP step to filter out baselines and channels
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DPPP step to filter out baselines and channels
/// @author Ger van Diepen

#ifndef DPPP_FILTER_H
#define DPPP_FILTER_H

#include "InputStep.h"

#include "../base/DPBuffer.h"
#include "../base/BaselineSelection.h"

namespace dp3 {
namespace common {
class ParameterSet;
}

namespace steps {
/// @brief DPPP step to filter out baselines and channels

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

class Filter : public Step {
 public:
  /// Default constructor.
  Filter();

  /// Construct the object for the given MS.
  /// Parameters are obtained from the parset using the given prefix.
  Filter(InputStep* input, const common::ParameterSet&, const string& prefix);

  /// Construct the object for the given MS and baseline selection.
  Filter(InputStep* input, const base::BaselineSelection&);

  virtual ~Filter();

  /// Process the next data chunk.
  /// When processed, it invokes the process function of the next step.
  virtual bool process(const base::DPBuffer&);

  /// Finish the processing of this step and subsequent steps.
  virtual void finish();

  /// Update the general info.
  virtual void updateInfo(const base::DPInfo&);

  /// Show the step parameters.
  virtual void show(std::ostream&) const;

  /// Show the timings.
  virtual void showTimings(std::ostream&, double duration) const;

  /// If needed, remove the deleted stations from the subtables
  /// and renumber the remaining stations.
  virtual void addToMS(const string& msName);

  /// Does the filter step has an actual selection?
  bool hasSelection() const { return itsDoSelect; }

  /// Get the indices of the selected baselines.
  const std::vector<unsigned int>& getIndicesBL() const { return itsSelBL; }

  /// Get the buffer.
  const base::DPBuffer& getBuffer() const { return itsBuf; }

 private:
  /// Create the mapping from old to new id (e.g. ANTENNA_ID).
  /// The removed ids get a mapping -1.
  casacore::Vector<int> createIdMap(
      common::rownr_t nrId,
      const casacore::Vector<common::rownr_t>& removedIds) const;

  /// Remove rows with deleted stations from a subtable.
  /// Renumber the ANTENNA_ID of the remaining rows.
  /// It fills nrId with the original number of rows in the subtable
  /// and returns the vector of removed row numbers.
  casacore::Vector<common::rownr_t> renumberSubTable(
      const casacore::Table& ms, const casacore::String& name,
      const casacore::String& colName,
      const casacore::Vector<common::rownr_t>& removedAnt,
      const casacore::Vector<int>& antMap, common::rownr_t& nrId) const;

  InputStep* itsInput;
  string itsName;
  base::DPBuffer itsBuf;
  base::DPBuffer itsBufTmp;
  casacore::String itsStartChanStr;  ///< startchan expression
  casacore::String itsNrChanStr;     ///< nchan expression
  bool itsRemoveAnt;                 ///< Remove from ANTENNA table?
  base::BaselineSelection itsBaselines;
  unsigned int itsStartChan;
  std::vector<unsigned int> itsSelBL;  ///< Index of baselines to select
  bool itsDoSelect;                    ///< Any selection?
  common::NSTimer itsTimer;
};

}  // namespace steps
}  // namespace dp3

#endif
