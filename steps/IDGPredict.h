// Copyright (C) 2023 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DP3_STEPS_IDGPREDICT_H_
#define DP3_STEPS_IDGPREDICT_H_

#include <complex>
#include <functional>
#include <string>
#include <vector>
#include <utility>

#ifdef HAVE_IDG
#include <idg-api.h>
#endif

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <EveryBeam/aterms/atermbase.h>

#include <schaapcommon/facets/facetimage.h>
#include <schaapcommon/facets/facet.h>

#include <aocommon/fits/fitsreader.h>
#include <aocommon/uvector.h>

#include <dp3/steps/Step.h>

#include "../common/ParameterSet.h"
#include "../common/Timer.h"

namespace dp3 {
namespace steps {

class IDGPredict : public ModelDataStep {
 public:
  IDGPredict(const common::ParameterSet& parset, const std::string& prefix,
             const std::vector<aocommon::FitsReader>& readers,
             std::vector<schaapcommon::facets::Facet>&& facets,
             const std::string& ds9_regions_file = "");

  IDGPredict(const common::ParameterSet&, const std::string& prefix);

  common::Fields getRequiredFields() const override { return kUvwField; }

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
  std::vector<aocommon::xt::UTensor<std::complex<float>, 3>> Predict(
      size_t direction);

  bool IsStarted() const;

  base::Direction GetFirstDirection() const override;

  void SetBufferSize(size_t nTimesteps);
  size_t GetBufferSize() const;

  /// Read the fits files (nterms) for the idg prediction.
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

#ifdef HAVE_IDG
 private:
  /// Reads the model image for each FitsReader in readers_.
  /// @return A vector with the model image data for each image.
  std::vector<aocommon::UVector<float>> GetModelImages();

  /// Initializes IDG buffersets for all directions and terms.
  void StartIDG();

  /// Initializes the aterms_ and aterm_values_ lists.
  void InitializeATerms(const casacore::MeasurementSet& ms);

  /// Calculates the ATerms IDG should use.
  /// @param direction Direction index.
  /// @return The ATerms for the given direction.
  aocommon::UVector<std::complex<float>> GetAtermValues(size_t direction) const;

  /// Creates a vector with uvw pointers. Raises max_w_ if needed.
  /// @return A vector with pointers to the uvw values of the input buffers.
  std::vector<const double*> InitializeUVWs();

  /// Initializes output buffers and fills them with IDG predictions.
  /// @param direction Direction index.
  /// @param uvws uvw pointers from InitializeUVWs.
  /// @param term_data Buffer for storing results of non-first terms.
  ///        The returned result buffers have the results for the first term.
  /// @return Result buffer.
  std::vector<aocommon::xt::UTensor<std::complex<float>, 3>>
  ComputeVisibilities(size_t direction, const std::vector<const double*>& uvws,
                      std::complex<float>* term_data) const;

  /// Applies phase shift and polynomial term corrections to computed
  /// visibilities.
  /// @param result Result buffer, as computed by ComputeVisibilities().
  /// @param term_data Buffer that has results for non-first terms.
  void CorrectVisibilities(
      std::vector<aocommon::xt::UTensor<std::complex<float>, 3>>& result,
      const std::complex<float>* term_data);

  /// Returns the amount of buffers that can be used by this step.
  /// If multiple IDG predicts have to run simultaneously, you can update the
  /// buffer size by using GetBufferSize and SetBufferSize respectively.
  size_t GetAllocatableBuffers(size_t memory);

  /// @return The number of items in a subgrid, for one antenna.
  size_t GetSubgridCount(size_t direction) const;

  std::string name_;
  common::ParameterSet parset_;

  std::vector<schaapcommon::facets::FacetImage> images_;
  std::vector<std::unique_ptr<idg::api::BufferSet>> buffersets_;

  struct FacetMetaData {
    FacetMetaData(double _dl, double _dm, double _dp)
        : dl(_dl), dm(_dm), dp(_dp) {}
    double dl, dm, dp;
  };
  std::vector<FacetMetaData> meta_data_;

  std::vector<std::unique_ptr<base::DPBuffer>> buffers_;

  double ref_frequency_;
  double pixel_size_x_, pixel_size_y_;
  std::vector<aocommon::FitsReader> readers_;

  size_t buffer_size_;  ///< Number of DPBuffers to keep before calling flush

  std::vector<std::size_t> ant1_;  ///< Contains only the used antennas
  std::vector<std::size_t> ant2_;  ///< Contains only the used antennas

  common::NSTimer timer_;
  double max_w_;
  double max_baseline_;
  std::vector<std::pair<double, double>> directions_;
  bool save_facets_;  ///< Write the facets?

  // Required for aterms (beam)
  std::vector<std::unique_ptr<everybeam::aterms::ATermBase>> aterms_;
  /// For every aterm, stores the values returned from aterm.Calculate(). This
  /// member is mutable since it acts as a cache: If aterm.Calculate() returns
  /// false, GetAtermValues() should use the old / cached values.
  mutable std::vector<aocommon::UVector<std::complex<float>>> aterm_values_;
#endif
};

}  // namespace steps
}  // namespace dp3

#endif  // header guard
