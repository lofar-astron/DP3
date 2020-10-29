// Copyright (C) 2020
// ASTRON (Netherlands Institute for Radio Astronomy)
// P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
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

#ifndef IDG_PREDICT_H
#define IDG_PREDICT_H

#ifdef HAVE_IDG
#include "FacetImage.h"

#include <idg-api.h>

#include <aocommon/fits/fitsreader.h>
#endif

#include "../Common/ParameterSet.h"
#include "../DPPP/DPStep.h"

#include <complex>
#include <functional>
#include <string>
#include <vector>
#include <utility>

using aocommon::FitsReader;

namespace DP3 {
namespace DPPP {

class IDGPredict : public DPStep {
 public:
  IDGPredict(
      DPInput* input, const ParameterSet&, const string& prefix,
      std::pair<std::vector<FitsReader>, std::vector<aocommon::UVector<double>>>
          readers,
      std::vector<Facet>& facets, const std::string& ds9_regions_file = "");

  IDGPredict(DPInput* input, const ParameterSet&, const string& prefix);

  void updateInfo(const DPInfo& info) override;

  /** Add a buffer to the IDG predictor, for use in Predict(), later. Calls
   * FlushBuffers if the buffer is full. */
  bool process(const DPBuffer& buffer) override;

  void finish() override;

  void show(std::ostream&) const override;

  /// Show the timings.
  void showTimings(std::ostream&, double duration) const override;

  /** Flush all available buffers to IDG. */
  void FlushBuffers();

  /**
   * Predict visibilities for added buffers in a given direction.
   * @param direction Index for the requested direction.
   * @return Buffers with the predicted visibilities. For each buffer added
   *         with process(), there is one corresponding output buffer.
   */
  std::vector<DPBuffer> Predict(size_t direction);

  bool IsStarted() const;

  const std::vector<std::pair<double, double>>& GetDirections() const;

  void SetBufferSize(size_t nTimesteps);
  const size_t GetBufferSize() const { return buffer_size_; }

  /// Read the tif files (nterms) for the idg prediction.
  static std::pair<std::vector<FitsReader>,
                   std::vector<aocommon::UVector<double>>>
  GetReaders(const std::vector<std::string>& fits_model_files);

  /// Get the facets from a region file and create the image models with the
  /// given image size
  static void GetFacets(std::vector<Facet>& facets_out,
                        const std::string& ds9_regions_file, const double ra,
                        const double dec, const double pixel_size_x,
                        const double pixel_size_y, const size_t full_width,
                        const size_t full_height);

  /// Get the facets from a region file and use readers to create the image
  /// models.
  static void GetFacets(std::vector<Facet>& facets_out,
                        const std::string& ds9_regions_file,
                        const FitsReader& reader);

 private:
  void StartIDG();
#ifdef HAVE_IDG

  std::vector<const double*> InitializeUVWs();

  std::vector<DPBuffer> ComputeVisibilities(
      size_t direction, const std::vector<const double*>& uvws,
      std::complex<float>* term_data);

  double ComputePhaseShiftFactor(const double* uvw, size_t direction) const;

  void CorrectVisibilities(const std::vector<const double*>& uvws,
                           std::vector<DPBuffer>& result,
                           const std::complex<float>* term_data,
                           size_t direction);

  /// Return the amount of buffers that can be used by this step.
  /// If multiple IDG predicts are ran simultaneously, you can update the
  /// buffer size by using GetBufferSize and SetBufferSize respectively.
  size_t GetAllocatableBuffers(size_t memory);

  std::string name_;

  std::vector<FacetImage> images_;
  std::vector<std::unique_ptr<idg::api::BufferSet>> buffersets_;
  struct FacetMetaData {
    FacetMetaData(double _dl, double _dm, double _dp)
        : dl(_dl), dm(_dm), dp(_dp) {}
    double dl, dm, dp;
  };
  std::vector<FacetMetaData> meta_data_;

  std::vector<DPBuffer> buffers_;

  double ref_frequency_;
  double pixel_size_x_, pixel_size_y_;
  std::vector<FitsReader> readers_;

  // Currently unused, will be useful when IDGPredict is a DPStep.
  size_t buffer_size_;

  DPInput* input_;
  std::vector<std::size_t> ant1_;  // Contains only the used antennas
  std::vector<std::size_t> ant2_;  // Contains only the used antennas

  NSTimer timer_;
  double max_w_;
  double max_baseline_;
  std::vector<std::pair<double, double>> directions_;
  bool save_facets_;  ///< Write the facets?
#endif
};

}  // namespace DPPP
}  // namespace DP3

#endif  // header guard
