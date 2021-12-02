// MSWriter.h: DPPP step writing to an MS
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DPPP step writing to an MS
/// @author Ger van Diepen

#ifndef DPPP_MSWRITER_H
#define DPPP_MSWRITER_H

#include "Step.h"

#include "../base/StManParsetKeys.h"

#include <casacore/tables/Tables/ColumnDesc.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/ArrayColumn.h>

namespace casacore {
class Table;
}

namespace dp3 {
namespace common {
class ParameterSet;
}

namespace steps {
class InputStep;

/// @brief DPPP step writing to an MS

/// This class is a Step creating a new MeasurementSet and writing
/// all data in it.
/// Most meta information (subtables and meta columns in main table) is
/// copied from the input MeasurementSet given by the MSReader object.
/// <br>
/// In principle the new MS uses the same storage managers as used in the
/// input MS, but in case of an MS stored with LofarStMan it will use the
/// optimal storage managers (ISM for slowly varying meta data, TSM for
/// bulk data, SSM for others).
///
/// The SPECTRAL_WINDOW table will be changed to reflect the channels
/// being used or averaged.
/// The OBSERVATION table will be updated for the correct start and end time.
/// The HISTORY table gets an entry containing the parset values and the
/// DPPP version.

class MSWriter : public Step {
 public:
  explicit MSWriter(InputStep* reader, const std::string& outName,
                    const common::ParameterSet&, const std::string& prefix);

  virtual ~MSWriter();

  /// Process the next data chunk.
  /// It returns false when at the end.
  virtual bool process(const base::DPBuffer&);

  /// Finish the processing of this step and subsequent steps.
  virtual void finish();

  /// Show the step parameters.
  virtual void show(std::ostream&) const;

  /// Add some data to the MeasurementSet written/updated.
  /// Calls addToMS from the previous step, with the current output msname.
  virtual void addToMS(const string&);

  /// Update the general info.
  virtual void updateInfo(const base::DPInfo&);

  /// Show the timings.
  virtual void showTimings(std::ostream&, double duration) const;

  /// Write the parset info into the HISTORY table of the MS.
  static void writeHistory(casacore::Table& ms,
                           const common::ParameterSet& parset);

  static void updateBeam(const std::string& outName,
                         const std::string& outColName,
                         const base::DPInfo& info);

 private:
  /// Create an array column description and add to table with given
  /// stoage manager (if given).
  void makeArrayColumn(casacore::ColumnDesc desc,
                       const casacore::IPosition& shape,
                       casacore::DataManager* dm, casacore::Table& table,
                       bool makeDirectColumn = false);

  /// Create the MS by cloning all subtables from the input MS.
  /// All output columns in the main table are using normal storage managers.
  /// The SPECTRAL_WINDOW table is adapted as needed.
  void createMS(const std::string& outName, const base::DPInfo& info,
                unsigned int tileSize, unsigned int tileNChan);

  /// Update the SPECTRAL_WINDOW table for averaged channels.
  void updateSpw(const string& outName, const base::DPInfo& info);

  /// Update the OBSERVATION table with the correct start and end time.
  void updateObs(const string& outName);

  /// Update the FIELD table with the new phase center.
  void updatePhaseCentre(const string& outName, const base::DPInfo& info);

  /// Write the data, flags, etc.
  void writeData(casacore::Table& out, const base::DPBuffer& buf);

  /// Write the full resolution flags (flags before any averaging).
  void writeFullResFlags(casacore::Table& out, const base::DPBuffer& buf);

  /// Write all meta data columns for a time slot (ANTENNA1, etc.)
  void writeMeta(casacore::Table& out, const base::DPBuffer& buf);

  /// Copy meta data columns for a time slot (ANTENNA1, etc.)
  /// It also copies all time info if possible.
  void copyMeta(const casacore::Table& in, casacore::Table& out,
                bool copyTimeInfo);

  /// Copy the contents of a scalar column.
  template <typename T>
  void fillSca(const T& value, casacore::Table& out,
               const casacore::String& columnName) {
    casacore::ScalarColumn<T> outCol(out, columnName);
    outCol.fillColumn(value);
  }

  /// Copy the contents of an array column.
  template <typename T>
  void fillArr(const casacore::Array<T>& value, casacore::Table& out,
               const casacore::String& columnName) {
    casacore::ArrayColumn<T> outCol(out, columnName);
    outCol.fillColumn(value);
  }

  /// Copy the contents of a scalar column.
  template <typename T>
  void copySca(const casacore::Table& in, casacore::Table& out,
               const casacore::String& columnName) {
    casacore::ROScalarColumn<T> inCol(in, columnName);
    casacore::ScalarColumn<T> outCol(out, columnName);
    outCol.putColumn(inCol.getColumn());
  }

  /// Copy the contents of an array column.
  template <typename T>
  void copyArr(const casacore::Table& in, casacore::Table& out,
               const casacore::String& columnName) {
    casacore::ROArrayColumn<T> inCol(in, columnName);
    casacore::ArrayColumn<T> outCol(out, columnName);
    outCol.putColumn(inCol.getColumn());
  }

  InputStep* itsReader;
  string itsName;
  string itsOutName;
  base::DPBuffer itsBuffer;
  casacore::Table itsMS;
  common::ParameterSet itsParset;  ///< parset for writing history
  casacore::String itsDataColName;
  casacore::String itsFlagColName;
  casacore::String itsWeightColName;
  double itsInterval;
  bool itsOverwrite;  ///< Overwrite an existing output MS?
  bool itsCopyCorrData;
  bool itsCopyModelData;
  bool itsWriteFullResFlags;
  unsigned int itsTileSize;
  unsigned int itsTileNChan;
  unsigned int itsNrCorr;
  unsigned int itsNrChan;
  unsigned int itsNrBl;
  unsigned int itsNrTimes;
  unsigned int itsNChanAvg;      ///< nr of channels in input averaged to 1
  unsigned int itsNTimeAvg;      ///< nr of times in input averaged to 1
  unsigned int itsNrTimesFlush;  ///< flush every N time slots (0=no flush)
  unsigned int itsNrDone;        ///< nr of time slots written
  std::string itsVdsDir;         ///< directory where to put VDS file
  std::string itsClusterDesc;    ///< name of clusterdesc file
  common::NSTimer itsTimer;
  base::StManParsetKeys itsStManKeys;
};

}  // namespace steps
}  // namespace dp3

#endif
