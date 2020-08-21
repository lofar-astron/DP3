// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
//
// This file is part of the LOFAR software suite.
// The LOFAR software suite is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The LOFAR software suite is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.

/// @file
/// @brief Step for writing BDA data to an MS.
/// @author Maik Nijhuis

#ifndef DPPP_MSBDAWRITER_H
#define DPPP_MSBDAWRITER_H

#include <casacore/tables/Tables/Table.h>
#include <map>

#include "DPStep.h"
#include "MSReader.h"
#include "../Common/ParameterSet.h"

using casacore::Int;
using casacore::Table;

namespace DP3 {
namespace DPPP {

class MSBDAWriter : public DPStep {
 public:
  MSBDAWriter(DPInput*, const string&, const ParameterSet&, const string&);

  ~MSBDAWriter() override;

  void updateInfo(const DPInfo&) override;

  virtual bool process(std::unique_ptr<BDABuffer>) override;

  void finish() override;

  void show(std::ostream&) const override;

 private:
  /// Create the MS by cloning all subtables from the input MS.
  /// All output columns in the main table are using normal storage managers.
  void CreateMS();

  void CreateMainTable();

  ///  Add the BDA_TIME_AXIS table to the measurement set
  /// if it does not already exist.
  void CreateBDATimeAxis();

  ///  Add the BDA_TIME_FACTOR table to the measurement set
  /// if it does not already exist.
  void CreateBDATimeFactor();

  /// Add the BDA_FREQ_AXIS_ID and BDA_SET_ID columns to SPECTRAL_WINDOW.
  void CreateMetaDataFrequencyColumns();

  /// Write a metadata row to BDA_TIME_AXIS
  /// metadata of the baselines to BDA_TIME_FACTOR
  /// and write to the metadata columns of SPECTRAL_WINDOW.
  /// If an entry already exists, nothing is written.
  void WriteMetaData();

  /// Write all the baselines to the BDA_TIME_FACTOR table.
  void WriteTimeFactorRows(const Int&, unsigned int&, unsigned int&);

  /// Write a row in the BDA_TIME_AXIS table.
  void WriteTimeAxisRow(const Int&, const unsigned int&, const unsigned int&);

  /// Overwrite the SPECTRAL_WINDOW and DATA_DESCRIPTION tables.
  void OverwriteSubTables(const Int&);

 private:
  DPInput* reader_;
  const std::string outName_;
  const ParameterSet parset_;
  const std::string prefix_;
  const bool overwrite_;

  unsigned int ncorr_;
  unsigned int nbl_;
  std::map<std::size_t, unsigned int> nchanToDescId;
  casacore::Table ms_;
};

}  // namespace DPPP
}  // namespace DP3

#endif
