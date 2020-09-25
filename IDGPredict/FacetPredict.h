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

#include <complex>
#include <functional>
#include <string>
#include <vector>

class FacetPredict {
 public:
  FacetPredict(const std::vector<std::string> fitsModelFiles,
               const std::string& ds9RegionsFile);

  void SetMSInfo(std::vector<std::vector<double>>&& bands, size_t nr_stations);

  bool IsStarted() const { return !_buffersets.empty(); }

  void StartIDG(bool saveFacets);

  void RequestPredict(size_t direction, size_t dataDescId, size_t rowId,
                      size_t timeIndex, size_t antenna1, size_t antenna2,
                      const double* uvw);

  std::function<void(size_t /*row*/, size_t /*direction*/,
                     size_t /*dataDescId*/,
                     const std::complex<float>* /*values*/)>
      PredictCallback;

  size_t NDirections() const { return _images.size(); }

  std::pair<double, double> Direction(size_t facet) const {
    return _directions[facet];
  }

  void Flush();

  void SetBufferSize(size_t nTimesteps) { _bufferSize = nTimesteps; }

 private:
  void computePredictionBuffer(size_t dataDescId, size_t direction);

  constexpr static double c() { return 299792458.0L; }

  std::vector<FacetImage> _images;
  std::vector<std::unique_ptr<idg::api::BufferSet>> _buffersets;
  struct FacetMetaData {
    double dl, dm, dp;
    bool isInitialized;
    size_t rowIdOffset;
    std::vector<double> uvws;
  };
  std::vector<FacetMetaData> _metaData;

  size_t _fullWidth, _fullHeight;
  double _refFrequency;
  double _pixelSizeX, _pixelSizeY;
  std::vector<FitsReader> _readers;
  double _padding;
  size_t _bufferSize;

  /// MS info
  double _maxW;
  std::vector<std::vector<double>> _bands;
  size_t _nr_stations;
  double _maxBaseline;
  std::vector<std::pair<double, double>> _directions;
};

#else  // HAVE_IDG

#include <complex>
#include <functional>
#include <string>
#include <vector>

class FacetPredict {
 public:
  FacetPredict(const std::vector<std::string>&, const std::string&) {
    notCompiled();
  }

  void SetMSInfo(std::vector<std::vector<double>>&& bands, size_t nr_stations) {
    notCompiled();
  }

  bool IsStarted() const {
    notCompiled();
    return false;
  }

  void StartIDG(bool) { notCompiled(); }

  void RequestPredict(size_t, size_t, size_t, size_t, size_t, size_t,
                      const double*) {
    notCompiled();
  }

  std::function<void(size_t, size_t, size_t, const std::complex<float>*)>
      PredictCallback;

  size_t NDirections() const {
    notCompiled();
    return 0;
  }

  std::pair<double, double> Direction(size_t) const {
    notCompiled();
    return std::pair<double, double>();
  }

  void Flush() { notCompiled(); }

  void SetBufferSize(size_t) { notCompiled(); }

 private:
  void notCompiled() const {
    throw std::runtime_error(
        "Facet prediction is not available, because DP3 was not compiled with "
        "IDG support");
  }
};

#endif  // HAVE_IDG

#endif  // header guard
