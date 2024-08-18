// MultiMSReader.h: DP3 step reading from multiple MSs
// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

/// @file
/// @brief DP3 step reading from multiple MSs
/// @author Ger van Diepen

#ifndef DP3_STEPS_MULTIMSREADER_H_
#define DP3_STEPS_MULTIMSREADER_H_

#include <casacore/tables/Tables/TableIter.h>
#include <casacore/tables/Tables/RefRows.h>
#include <casacore/casa/Arrays/Slicer.h>

#include "../base/UVWCalculator.h"
#include "../base/FlagCounter.h"
#include "../common/ParameterSet.h"

#include "MSReader.h"
#include "ResultStep.h"

namespace dp3 {
namespace steps {
/// @brief DP3 step reading from multiple MSs

/// This class is an InputStep step reading the data from multiple
/// MeasurementSets (MSs) which can be used when the total number of channels
/// has been distributed over multiple MSs. It is therefore important that the
/// shape of all data inside the MSs (visibilities, flags, weights, etc.) is
/// identical except for the channel dimension.
///
/// Similar to the MSReader step, the object is constructed from the
/// 'msin' keywords in the parset file. These keywords are identical to those
/// in the MSReader because this class creates a vector of several MSReader
/// steps, where each MSReader is responsible for reading one of the MSs.
/// So, refer to the documentation of the MSReader for the definitions of the
/// 'msin' keywords that can be used for this class.
/// In addition to those keywords, the following can be given:
/// <ul>
///  <li> msin.orderms: order the MSs on frequency? If yes, all MSs must exist,
///           otherwise they cannot be ordered. If no, the MSs must be given
///           in order of frequency [yes]
///  <li> msin.missingdata: allow a non-existing data column in an MS? [no]
/// </ul>

class MultiMSReader final : public MSReader {
 public:
  /// Construct the object for the given MS.
  /// Parameters are obtained from the parset using the given prefix.
  MultiMSReader(const std::vector<string>& msNames,
                const common::ParameterSet& parset, const std::string& prefix);

  ~MultiMSReader() override;

  /// Process the next data chunk.
  /// It returns false when at the end.
  bool process(std::unique_ptr<base::DPBuffer> buffer) override;

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

  /// Set which fields must be read.
  void setFieldsToRead(const dp3::common::Fields& fields) override;

  /// Get the name of the first MS.
  std::string msName() const override;

 private:
  /// Validate that all bands have matching properties.
  /// @throw std::runtime_error If a band has non-matching properties.
  void ValidateBands();

  /// Handle the band info when all MSs are present.
  void HandleBands();

  /// Sort the bands (MSs) in order of frequency.
  void SortBands();

  /// Fill the band info where some MSs are missing.
  void FillBands();

  /// Reads the weights into 'buffer'
  void getWeights(std::unique_ptr<base::DPBuffer>& buffer);

  bool itsOrderMS;  ///< sort multi MS in order of freq?
  int itsFirst;     ///< first valid MSReader (<0 = none)
  int itsNMissing;  ///< nr of missing MSs
  struct Reader {
    std::string name;
    std::shared_ptr<MSReader> ms_reader;
    std::shared_ptr<ResultStep> result;
  };
  std::vector<Reader> readers_;
  unsigned int itsFillNChan;  ///< nr of channels for missing MSs
};

}  // namespace steps
}  // namespace dp3

#endif
