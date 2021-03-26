// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief Step for writing BDA data to an MS.
/// @author Maik Nijhuis

#ifndef DPPP_MSBDAWRITER_H
#define DPPP_MSBDAWRITER_H

#include <casacore/tables/Tables/Table.h>
#include <map>

#include "Step.h"
#include "MSReader.h"
#include "../common/ParameterSet.h"

namespace dp3 {
namespace steps {

class MSBDAWriter : public Step {
 public:
  MSBDAWriter(InputStep*, const std::string&, const common::ParameterSet&,
              const std::string&);

  ~MSBDAWriter() override;

  void updateInfo(const base::DPInfo&) override;

  virtual bool process(std::unique_ptr<base::BDABuffer>) override;

  void finish() override;

  /// Add some data to the MeasurementSet written/updated.
  /// Calls addToMS from the previous step, with the current output msname.
  void addToMS(const std::string&) override;

  void show(std::ostream&) const override;

 private:
  /// Create the MS by cloning all subtables from the input MS.
  /// All output columns in the main table are using normal storage managers.
  void CreateMS();

  void CreateMainTable();

  ///  Add the BDA_TIME_AXIS table to the measurement set
  /// if it does not already exist.
  void CreateBDATimeAxis();

  ///  Add the BDA_FACTORS table to the measurement set
  /// if it does not already exist.
  void CreateBDATimeFactor();

  /// Add the BDA_FREQ_AXIS_ID and BDA_SET_ID columns to SPECTRAL_WINDOW.
  void CreateMetaDataFrequencyColumns();

  /// Write a metadata row to BDA_TIME_AXIS
  /// metadata of the baselines to BDA_FACTORS
  /// and write to the metadata columns of SPECTRAL_WINDOW.
  /// If an entry already exists, nothing is written.
  void WriteMetaData();

  /// Write all the baselines to the BDA_FACTORS table,
  /// while updating min_factor_time and max_factor_time.
  void WriteTimeFactorRows(casacore::Int bda_set_id,
                           unsigned int& min_factor_time,
                           unsigned int& max_factor_time);

  /// Write a row in the BDA_TIME_AXIS table.
  void WriteTimeAxisRow(casacore::Int bda_set_id, unsigned int min_factor_time,
                        unsigned int max_factor_time);

  /// Overwrite the SPECTRAL_WINDOW and DATA_DESCRIPTION tables.
  void OverwriteSubTables(casacore::Int bda_set_id);

 private:
  InputStep* reader_;
  const std::string outName_;
  const common::ParameterSet parset_;
  const std::string prefix_;
  const bool overwrite_;

  unsigned int ncorr_;
  unsigned int nbl_;
  std::map<std::size_t, unsigned int> nchanToDescId;
  casacore::Table ms_;
};

}  // namespace steps
}  // namespace dp3

#endif
