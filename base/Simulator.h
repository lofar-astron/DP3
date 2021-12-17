// Simulator.h: Compute visibilities for different model components types
// (implementation of ModelComponentVisitor).
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DPPP_SIMULATOR_H
#define DPPP_SIMULATOR_H

#include "Baseline.h"
#include "ModelComponent.h"
#include "ModelComponentVisitor.h"
#include "Direction.h"

#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/Cube.h>

namespace dp3 {
namespace base {

/// @{

typedef std::complex<double> dcomplex;

/// @brief Simulator for different kinds of model components
class Simulator : public ModelComponentVisitor {
 public:
  /**
   * @brief Construct a new Simulator object
   *
   * @param reference An existing BDABuffer.
   * @param baselines Vector of Baselines
   * @param freq Channel frequencies (Hz)
   * @param chanWidths Channel widths (Hz)
   * @param stationUVW Station UVW coordinates
   * @param buffer Output buffer, should be of shape (nCor, nFreq, nBaselines),
   * where nCor should be 1 if stokesIOnly is true, else 4
   * @param correctFreqSmearing Correct for frequency smearing
   * @param stokesIOnly Stokes I only, to avoid a loop over correlations
   */
  Simulator(const Direction& reference, size_t nStation,
            const std::vector<Baseline>& baselines,
            const casacore::Vector<double>& freq,
            const casacore::Vector<double>& chanWidths,
            const casacore::Matrix<double>& stationUVW,
            casacore::Cube<dcomplex>& buffer, bool correctFreqSmearing,
            bool stokesIOnly);

  // Note DuoMatrix is actually two T matrices
  // T: floating point type, ideally float, double, or long double.
  template <typename T>
  class DuoMatrix {
   public:
    DuoMatrix() : itsNRows(0){};

    DuoMatrix(size_t nrows, size_t ncols) { resize(nrows, ncols); }

    void resize(size_t nrows, size_t ncols) {
      itsNRows = nrows;
      itsData_real.resize(nrows * ncols);
      itsData_imag.resize(nrows * ncols);
    }

    size_t nRows() { return itsNRows; }
    size_t nCols() { return itsData_real.size() / itsNRows; }

    T& real(size_t row, size_t col) {
      return itsData_real[col * itsNRows + row];
    }
    T& imag(size_t row, size_t col) {
      return itsData_imag[col * itsNRows + row];
    }
    T* realdata() { return &itsData_real[0]; }
    T* imagdata() { return &itsData_imag[0]; }

   private:
    // Use separate storage for real/imag parts
    std::vector<T> itsData_real;
    std::vector<T> itsData_imag;
    size_t itsNRows;
  };

  void simulate(const ModelComponent::ConstPtr& component);

 private:
  virtual void visit(const PointSource& component);
  virtual void visit(const GaussianSource& component);

 private:
  Direction itsReference;
  size_t itsNStation, itsNBaseline, itsNChannel;
  bool itsCorrectFreqSmearing;
  bool itsStokesIOnly;
  const std::vector<Baseline> itsBaselines;
  const casacore::Vector<double> itsFreq;
  const casacore::Vector<double> itsChanWidths;
  const casacore::Matrix<double> itsStationUVW;
  casacore::Cube<dcomplex> itsBuffer;
  std::vector<double> itsStationPhases;
  DuoMatrix<double> itsShiftBuffer;
  DuoMatrix<double> itsSpectrumBuffer;
};

/// @}

}  // namespace base
}  // namespace dp3

#endif
