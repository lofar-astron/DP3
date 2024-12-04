// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_STEPS_MSBDAWRITER_H_
#define DP3_STEPS_MSBDAWRITER_H_

#include <map>

#include <casacore/tables/Tables/Table.h>

#include "../common/ParameterSet.h"

#include "OutputStep.h"

namespace dp3 {
namespace steps {

/**
 * @brief Step for writing BDA data to an MS.
 */
class MSBDAWriter : public OutputStep {
 public:
  explicit MSBDAWriter(const std::string& out_name,
                       const common::ParameterSet& parset,
                       const std::string& prefix);

  common::Fields getRequiredFields() const override {
    return kDataField | kFlagsField | kWeightsField | kUvwField;
  }

  void updateInfo(const base::DPInfo&) override;

  bool process(std::unique_ptr<base::BdaBuffer>) override;

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
  const std::string out_name_;
  const common::ParameterSet parset_;
  const std::string prefix_;
  const bool overwrite_;

  std::map<std::size_t, unsigned int> nchanToDescId;
  casacore::Table ms_;
};

}  // namespace steps
}  // namespace dp3

#endif
