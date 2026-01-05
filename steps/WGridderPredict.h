// Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_STEPS_WGRIDDERPREDICT_H_
#define DP3_STEPS_WGRIDDERPREDICT_H_

#include <complex>
#include <functional>
#include <string>
#include <vector>
#include <utility>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <schaapcommon/facets/facetimage.h>
#include <schaapcommon/facets/facet.h>

#include <aocommon/fits/fitsreader.h>
#include <aocommon/uvector.h>

#include "steps/Step.h"

#include "common/ParameterSet.h"
#include "common/Timer.h"

namespace dp3 {
namespace steps {

class WGridderPredict : public ModelDataStep {
 public:
  WGridderPredict(const common::ParameterSet& parset, const std::string& prefix,
                  std::vector<aocommon::FitsReader>&& readers,
                  std::vector<schaapcommon::facets::Facet>&& facets,
                  const std::string& ds9_regions_file = "");

  WGridderPredict(const common::ParameterSet&, const std::string& prefix);

  common::Fields getRequiredFields() const override { return kUvwField; }

  common::Fields getProvidedFields() const override {
    common::Fields fields;
    if (sum_facets_) {
      // The predicted visibilities go to the main data buffer
      fields |= kDataField;
    } else {
      // TODO(AST-1241): Handle these dependencies using Fields.
    }
    return fields;
  }

  void updateInfo(const base::DPInfo& info) override;

  /// Add a buffer to the IDG predictor, for use in Predict(), later. Calls
  /// flush if the buffer is full.
  bool process(std::unique_ptr<base::DPBuffer> buffer) override;

  void finish() override;

  void show(std::ostream&) const override;

  void showTimings(std::ostream&, double duration) const override;

  /// Process the data in all internal buffers using IDG, and send the results
  /// to the next step using its process() function.
  void flush();

  /// Predict visibilities for added buffers in a given direction.
  /// @param direction Index for the requested direction.
  /// @return Buffers with the predicted visibilities. For each buffer added
  ///         with process(), there is one corresponding output buffer.
  void Predict(size_t direction,
               std::vector<base::DPBuffer::DataType*>& destinations);

  base::Direction GetFirstDirection() const override;

  void SetBufferSize(size_t nTimesteps);
  size_t GetBufferSize() const;

  /// Read the fits files (nterms) for image based prediction.
  static std::vector<aocommon::FitsReader> GetReaders(
      const std::vector<std::string>& fits_model_files);

  /// Get the facets from a region file and create the image models with the
  /// given image size
  static std::vector<schaapcommon::facets::Facet> GetFacets(
      const std::string& ds9_regions_file, const double ra, const double dec,
      const double pixel_size_x, const double pixel_size_y,
      const size_t full_width, const size_t full_height);

  /// Get the facets from a region file and use readers to create the image
  /// models.
  static std::vector<schaapcommon::facets::Facet> GetFacets(
      const std::string& ds9_regions_file, const aocommon::FitsReader& reader);

 private:
  /// Returns the amount of buffers that can be used by this step.
  /// If multiple predicts have to run simultaneously, you can update the
  /// buffer size by using GetBufferSize and SetBufferSize respectively.
  size_t GetAllocatableBuffers(size_t memory);

  /// Create a vector of padded images from the readers
  std::vector<aocommon::Image> GetModelImages();

  std::string name_;
  common::ParameterSet parset_;

  std::vector<schaapcommon::facets::FacetImage> images_;

  struct FacetMetaData {
    FacetMetaData(double _dl, double _dm, double _dp)
        : dl(_dl), dm(_dm), dp(_dp) {}
    double dl, dm, dp;
  };
  std::vector<FacetMetaData> meta_data_;

  std::vector<std::unique_ptr<base::DPBuffer>> buffers_;

  double reference_frequency_;
  double pixel_size_x_, pixel_size_y_;
  std::vector<aocommon::FitsReader> readers_;

  size_t buffer_size_;  ///< Number of DPBuffers to keep before calling flush

  common::NSTimer timer_;
  std::vector<dp3::base::Direction> directions_;
  std::vector<std::string> direction_labels_;
  bool save_facets_;
  bool sum_facets_;
};

}  // namespace steps
}  // namespace dp3

#endif  // header guard
