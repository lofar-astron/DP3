// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DPPP_STMANPARSETKEYS_H
#define DPPP_STMANPARSETKEYS_H

#include "../common/ParameterSet.h"

#include <string>

#include <casacore/casa/Containers/Record.h>

#include <boost/algorithm/string/case_conv.hpp>

namespace dp3 {
namespace base {
struct StManParsetKeys {
  std::string storage_manager_name;
  /// Bits per data float, or 0 if data column is not compressed
  unsigned int dysco_data_bit_rate;
  /// Bits per weight float, or 0 if weight column is not compressed
  unsigned int dysco_weight_bit_rate;
  /// Distribution assumed for compression; e.g. "Uniform" or
  /// "TruncatedGaussian"
  std::string dysco_distribution;
  /// For truncated distributions, the truncation point (e.g. 3 for 3 sigma in
  /// TruncGaus).
  double dysco_dist_truncation;
  /// Kind of normalization; "AF", "RF" or "Row".
  std::string dysco_normalization;
  /// Deflate compression level for Sisco (only used for writing).
  int sisco_deflate_level;
  /// Prediction order: -1 for no prediction, 0 for compressing differences, 1
  /// for using linear prediction using 2 previous values, 2 for using quadratic
  /// prediction using 3 previous values.
  int sisco_predict_level;

  explicit StManParsetKeys(const common::ParameterSet& parset,
                           const std::string& prefix) {
    storage_manager_name = boost::to_lower_copy(parset.getString(
        prefix + "storagemanager",
        parset.getString(prefix + "storagemanager.name", std::string())));
    if (storage_manager_name == "dysco") {
      dysco_data_bit_rate =
          parset.getInt(prefix + "storagemanager.databitrate", 10);
      dysco_weight_bit_rate =
          parset.getInt(prefix + "storagemanager.weightbitrate", 12);
      dysco_distribution = parset.getString(
          prefix + "storagemanager.distribution", "TruncatedGaussian");
      dysco_dist_truncation =
          parset.getDouble(prefix + "storagemanager.disttruncation", 2.5);
      dysco_normalization =
          parset.getString(prefix + "storagemanager.normalization", "AF");
    }
    if (storage_manager_name == "sisco") {
      sisco_deflate_level =
          parset.getInt(prefix + "storagemanager.deflate_level", 9);
      sisco_predict_level =
          parset.getInt(prefix + "storagemanager.predict_level", 2);
    }
  }

  std::string_view GetStorageManagerClassName() const {
    if (storage_manager_name == "stokes_i")
      return "StokesIStMan";
    else if (storage_manager_name == "dysco")
      return "DyscoStMan";
    else if (storage_manager_name == "sisco")
      return "SiscoStMan";
    else if (!storage_manager_name.empty())
      throw std::runtime_error("Unknown storage manager specified: " +
                               storage_manager_name);
    else
      return {};
  }

  casacore::Record GetSpecification() const {
    if (storage_manager_name == "dysco")
      return GetDyscoSpec();
    else if (storage_manager_name == "sisco")
      return GetSiscoSpec();
    else
      return casacore::Record();
  }

  casacore::Record GetSiscoSpec() const {
    casacore::Record result;
    result.define("deflate_level", sisco_deflate_level);
    result.define("predict_level", sisco_predict_level);
    return result;
  }

  casacore::Record GetDyscoSpec() const {
    casacore::Record dyscoSpec;
    dyscoSpec.define("distribution", dysco_distribution);
    dyscoSpec.define("normalization", dysco_normalization);
    dyscoSpec.define("distributionTruncation", dysco_dist_truncation);
    /// DPPP uses bitrate of 0 to disable compression of the data/weight column.
    /// However, Dysco does not allow the data or weight bitrates to be set to
    /// 0, so we set the values to something different. The values are not
    /// actually used.
    unsigned int data_bit_rate = dysco_data_bit_rate;
    if (data_bit_rate == 0) data_bit_rate = 16;
    dyscoSpec.define("dataBitCount", data_bit_rate);
    unsigned int weight_bit_rate = dysco_weight_bit_rate;
    if (weight_bit_rate == 0) weight_bit_rate = 16;
    dyscoSpec.define("weightBitCount", weight_bit_rate);
    return dyscoSpec;
  }
};
}  // namespace base
}  // namespace dp3
#endif
