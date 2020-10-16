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

#ifndef FACET_PREDICT_H
#define FACET_PREDICT_H

#ifdef HAVE_IDG
#include "FacetImage.h"

#include <idg-api.h>

#include "FitsReader.h"
#endif

#include "../DPPP/DPStep.h"

#include <complex>
#include <functional>
#include <string>
#include <vector>

namespace DP3 {
namespace DPPP {

class FacetPredict {
 public:
  FacetPredict(DPInput& input, const std::vector<std::string>& fitsModelFiles,
               const std::string& ds9RegionsFile);

  void updateInfo(const DPInfo& info);

  /** Add a buffer to the facet predictor, for use in Predict(), later. */
  void AddBuffer(const DPBuffer& buffer) { buffers_.emplace_back(buffer); }

  /** Remove all previously added buffers. */
  void FlushBuffers() { buffers_.clear(); }

  /**
   * Predict visibilities for added buffers in a given direction.
   * @param direction Index for the requested direction.
   * @return Buffers with the predicted visibilities. For each buffer added
   *         with AddBuffer(), there is one corresponding output buffer.
   */
  std::vector<DPBuffer> Predict(size_t direction);

  bool IsStarted() const;

  void StartIDG(bool saveFacets);

  const std::vector<std::pair<double, double>>& GetDirections() const;

  void SetBufferSize(size_t nTimesteps);

#ifdef HAVE_IDG
 private:
  std::vector<const double*> InitializeUVWs();

  std::vector<DPBuffer> ComputeVisibilities(
      size_t direction, const std::vector<const double*>& uvws,
      std::complex<float>* term_data) const;

  double ComputePhaseShiftFactor(const double* uvw, size_t direction) const;

  void CorrectVisibilities(const std::vector<const double*>& uvws,
                           std::vector<DPBuffer>& result,
                           const std::complex<float>* term_data,
                           size_t direction) const;

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

  // Currently unused, will be useful when FacetPredict is a DPStep.
  size_t buffer_size_;

  DPInput& input_;
  DPInfo info_;
  std::vector<std::size_t> ant1_;  // Contains only the used antennas
  std::vector<std::size_t> ant2_;  // Contains only the used antennas

  NSTimer timer_;
  double max_w_;
  double max_baseline_;
  std::vector<std::pair<double, double>> directions_;
#endif
};

}  // namespace DPPP
}  // namespace DP3

#endif  // header guard
