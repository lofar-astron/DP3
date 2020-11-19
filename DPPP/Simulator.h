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
#include "Position.h"

#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/Cube.h>

namespace DP3 {
namespace DPPP {

/// @{

typedef std::complex<double> dcomplex;

/// @brief Compute visibilities for different model components types
/// (implementation of ModelComponentVisitor).
class Simulator : public ModelComponentVisitor {
 public:
  Simulator(const Position& reference, size_t nStation, size_t nBaseline,
            size_t nChannel, const casacore::Vector<Baseline>& baselines,
            const casacore::Vector<double>& freq,
            const casacore::Matrix<double>& uvw,
            casacore::Cube<dcomplex>& buffer, bool stokesIOnly = false);

  template <typename T>
  class Matrix {
   public:
    Matrix() : itsNRows(0){};

    Matrix(size_t nrows, size_t ncols) { resize(nrows, ncols); }

    void resize(size_t nrows, size_t ncols) {
      itsNRows = nrows;
      itsData.resize(nrows * ncols);
    }

    T& operator()(size_t row, size_t col) {
      return itsData[col * itsNRows + row];
    }

    T* data() { return &itsData[0]; }

   private:
    std::vector<T> itsData;
    size_t itsNRows;
  };

  void simulate(const ModelComponent::ConstPtr& component);

 private:
  virtual void visit(const PointSource& component);
  virtual void visit(const GaussianSource& component);

 private:
  Position itsReference;
  size_t itsNStation, itsNBaseline, itsNChannel;
  bool itsStokesIOnly;
  const casacore::Vector<Baseline> itsBaselines;
  const casacore::Vector<double> itsFreq;
  const casacore::Matrix<double> itsUVW;
  casacore::Cube<dcomplex> itsBuffer;
  Matrix<dcomplex> itsShiftBuffer;
  Matrix<dcomplex> itsSpectrumBuffer;
};

/// @}

}  // namespace DPPP
}  // namespace DP3

#endif
