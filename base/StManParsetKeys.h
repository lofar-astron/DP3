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
  std::string stManName;
  unsigned int dyscoDataBitRate;  ///< Bits per data float, or 0 if data column
                                  ///< is not compressed
  unsigned int dyscoWeightBitRate;  ///< Bits per weight float, or 0 if weight
                                    ///< column is not compressed
  std::string dyscoDistribution;    ///< Distribution assumed for compression;
                                    ///< e.g. "Uniform" or "TruncatedGaussian"
  double dyscoDistTruncation;  ///< For truncated distributions, the truncation
                               ///< point (e.g. 3 for 3 sigma in TruncGaus).
  std::string
      dyscoNormalization;  ///< Kind of normalization; "AF", "RF" or "Row".

  explicit StManParsetKeys(const common::ParameterSet& parset,
                           const std::string& prefix) {
    stManName = boost::to_lower_copy(parset.getString(
        prefix + "storagemanager",
        parset.getString(prefix + "storagemanager.name", std::string())));
    if (stManName == "dysco") {
      dyscoDataBitRate =
          parset.getInt(prefix + "storagemanager.databitrate", 10);
      dyscoWeightBitRate =
          parset.getInt(prefix + "storagemanager.weightbitrate", 12);
      dyscoDistribution = parset.getString(
          prefix + "storagemanager.distribution", "TruncatedGaussian");
      dyscoDistTruncation =
          parset.getDouble(prefix + "storagemanager.disttruncation", 2.5);
      dyscoNormalization =
          parset.getString(prefix + "storagemanager.normalization", "AF");
    }
  }

  casacore::Record GetDyscoSpec() const {
    casacore::Record dyscoSpec;
    dyscoSpec.define("distribution", dyscoDistribution);
    dyscoSpec.define("normalization", dyscoNormalization);
    dyscoSpec.define("distributionTruncation", dyscoDistTruncation);
    /// DPPP uses bitrate of 0 to disable compression of the data/weight column.
    /// However, Dysco does not allow the data or weight bitrates to be set to
    /// 0, so we set the values to something different. The values are not
    /// actually used.
    unsigned int dataBitRate = dyscoDataBitRate;
    if (dataBitRate == 0) dataBitRate = 16;
    dyscoSpec.define("dataBitCount", dataBitRate);
    unsigned int weightBitRate = dyscoWeightBitRate;
    if (weightBitRate == 0) weightBitRate = 16;
    dyscoSpec.define("weightBitCount", weightBitRate);
    return dyscoSpec;
  }
};
}  // namespace base
}  // namespace dp3
#endif
