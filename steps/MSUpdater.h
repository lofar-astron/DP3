// MSUpdater.h: DPPP step writing to an MS
// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_STEPS_MSUPDATER_H_
#define DP3_STEPS_MSUPDATER_H_

#include <casacore/tables/Tables/ColumnDesc.h>
#include <casacore/tables/Tables/RefRows.h>
#include <casacore/tables/Tables/Table.h>

#include "../base/StManParsetKeys.h"
#include "../common/ParameterSet.h"
#include "../common/Timer.h"

#include "OutputStep.h"

namespace dp3 {
namespace steps {

/// @brief DP3 step writing to an MS

/// This class updates the flags in an existing MeasurementSet.
/// Hardly anything is done in this class.
/// It uses function putFlags in MsReader to do the actual write.
//
/// Like MSWriter it adds an entry to the HISTORY table of the MS
/// containing the parset values and DPPP version.

class MSUpdater : public OutputStep {
 public:
  explicit MSUpdater(std::string msName, const common::ParameterSet& parset,
                     const std::string& prefix, bool writeHistory = true);

  common::Fields getRequiredFields() const override;

  void SetFieldsToWrite(const common::Fields& fields) override;

  /// Process the next data chunk.
  /// It returns false when at the end.
  bool process(std::unique_ptr<base::DPBuffer> buffer) override;

  /// Finish the processing of this step and subsequent steps.
  void finish() override;

  /// Update the general info.
  void updateInfo(const base::DPInfo&) override;

  /// Show the step parameters.
  void show(std::ostream&) const override;

  /// Show the timings.
  void showTimings(std::ostream&, double duration) const override;

  /// Tests if an update of the buffer described in info to the MS msName
  /// is possible. When throwError is true, it will throw an error with a
  /// descriptive string before returning false
  static bool updateAllowed(const base::DPInfo& info, casacore::String msName,
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
  bool addColumn(const std::string& colname, const casacore::DataType dataType,
                 const casacore::ColumnDesc& cd);

  std::string itsName;
  std::string itsMSName;
  casacore::Table itsMS;
  const common::ParameterSet& itsParset;
  base::DPBuffer itsBuffer;
  std::string itsDataColName;
  std::string itsFlagColName;
  std::string itsWeightColName;
  unsigned int itsNrTimesFlush;  ///< flush every N time slots (0=no flush)
  unsigned int itsNrDone;        ///< nr of time slots written
  bool itsDataColAdded;          ///< has data column been added?
  bool itsFlagColAdded;          ///< has weight column been added?
  bool itsWeightColAdded;        ///< has weight column been added?
  bool itsWriteHistory;          ///< Should history be written?
  common::NSTimer itsTimer;
  unsigned int itsTileSize;
  base::StManParsetKeys itsStManKeys;
};

}  // namespace steps
}  // namespace dp3

#endif
