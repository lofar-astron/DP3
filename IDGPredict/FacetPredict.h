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

#include "../DPPP/DPInfo.h"

#include <complex>
#include <functional>
#include <string>
#include <vector>

class FacetPredict {
 public:
  using PredictCallback = std::function<void(
      size_t /*row*/, size_t /*direction*/, size_t /*data_desc_id*/,
      const std::complex<float>* /*values*/)>;

  FacetPredict(const std::vector<std::string> fitsModelFiles,
               const std::string& ds9RegionsFile, PredictCallback&& callback);

  void updateInfo(const DP3::DPPP::DPInfo& info);

  bool IsStarted() const;

  void StartIDG(bool saveFacets);

  void RequestPredict(size_t direction, size_t rowId, size_t timeIndex,
                      size_t antenna1, size_t antenna2, const double* uvw);

  const std::vector<std::pair<double, double>>& GetDirections() const;

  void Flush(std::size_t direction);

  void SetBufferSize(size_t nTimesteps);

#ifdef HAVE_IDG
 private:
  void computePredictionBuffer(size_t direction);

  constexpr static double c() { return 299792458.0L; }

  constexpr double wavelength(size_t channel) const {
    return c() / info_.chanFreqs()[channel];
  }

  PredictCallback predict_callback_;
  std::vector<FacetImage> images_;
  std::vector<std::unique_ptr<idg::api::BufferSet>> buffersets_;
  struct FacetMetaData {
    double dl;
    double dm;
    double dp;
    bool is_initialized;
    size_t row_id_offset;
    std::vector<double> uvws;
  };
  std::vector<FacetMetaData> meta_data_;

  std::size_t full_width_;
  std::size_t full_height_;
  double ref_frequency_;
  double pixel_size_x_;
  double pixel_size_y_;
  std::vector<FitsReader> readers_;
  double padding_;
  size_t buffer_size_;

  /// MS info
  DP3::DPPP::DPInfo info_;
  double max_w_;
  double max_baseline_;
  std::vector<std::pair<double, double>> directions_;
#endif  // HAVE_IDG
};

#endif  // FACET_PREDICT_H
