// MSUpdater.h: DPPP step writing to an MS
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DPPP step writing to an MS
/// @author Ger van Diepen

#ifndef DPPP_MSUPDATER_H
#define DPPP_MSUPDATER_H

#include "DPStep.h"
#include "StManParsetKeys.h"

#include <casacore/tables/Tables/Table.h>

namespace casacore {
class ColumnDesc;
class RefRows;
}  // namespace casacore

namespace DP3 {

class ParameterSet;

namespace DPPP {
class DPInput;

/// @brief DPPP step writing to an MS

/// This class updates the flags in an existing MeasurementSet.
/// Hardly anything is done in this class.
/// It uses function putFlags in MSReader to do the actual write.
//
/// Like MSWriter it adds an entry to the HISTORY table of the MS
/// containing the parset values and DPPP version.

class MSUpdater : public DPStep {
 public:
  MSUpdater(DPInput* reader, casacore::String msName,
            const ParameterSet& parset, const std::string& prefix,
            bool writeHistory = true);

  virtual ~MSUpdater();

  /// Process the next data chunk.
  /// It returns false when at the end.
  virtual bool process(const DPBuffer&);

  /// Finish the processing of this step and subsequent steps.
  virtual void finish();

  /// Update the general info.
  virtual void updateInfo(const DPInfo&);

  /// Add some data to the MeasurementSet written/updated.
  /// Calls addToMS from the previous step, with the current output msname.
  virtual void addToMS(const string&);

  /// Show the step parameters.
  virtual void show(std::ostream&) const;

  /// Show the timings.
  virtual void showTimings(std::ostream&, double duration) const;

  /// Tests if an update of the buffer described in info to the MS msName
  /// is possible. When throwError is true, it will throw an error with a
  /// descriptive string before returning false
  static bool updateAllowed(const DPInfo& info, casacore::String msName,
                            bool throwError = true);

 private:
  /// Write the flags at the given row numbers.
  void putFlags(const casacore::RefRows& rowNrs,
                const casacore::Cube<bool>& flags);

  /// Write the weights at the given row numbers
  void putWeights(const casacore::RefRows& rowNrs,
                  const casacore::Cube<float>& weights);

  /// Write the data at the given row numbers.
  void putData(const casacore::RefRows& rowNrs,
               const casacore::Cube<casacore::Complex>& data);

  /// If not existing yet, add the column specified by colname.
  /// Column will containt arrays of type datatype.
  /// If the column has been added, this function returns true
  bool addColumn(const string& colname, const casacore::DataType dataType,
                 const casacore::ColumnDesc& cd);

  DPInput* itsReader;
  string itsName;
  casacore::String itsMSName;
  casacore::Table itsMS;
  const ParameterSet& itsParset;
  DPBuffer itsBuffer;
  casacore::String itsDataColName;
  casacore::String itsWeightColName;
  unsigned int itsNrTimesFlush;  ///< flush every N time slots (0=no flush)
  bool itsWriteData;
  bool itsWriteWeights;
  bool itsWriteFlags;
  unsigned int itsNrDone;  ///< nr of time slots written
  bool itsDataColAdded;    ///< has data column been added?
  bool itsWeightColAdded;  ///< has weight column been added?
  bool itsWriteHistory;    ///< Should history be written?
  NSTimer itsTimer;
  unsigned int itsTileSize;
  StManParsetKeys itsStManKeys;
};

}  // namespace DPPP
}  // namespace DP3

#endif
