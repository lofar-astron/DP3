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
#include <dp3/base/Direction.h>

#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/Cube.h>

namespace dp3 {
namespace base {

/// @{

typedef std::complex<double> dcomplex;

/**
 * @brief Simulator to compute visibilities given a sky model
 *
 * This class computes visibilities given model components in the sky.
 * Effectively, it evaluates the equation:
 *
 * \f[ V(u, v, w) =
 * \iint I(\ell, m) e^{2\pi i(u,v,w)\cdot(\ell,m,n)}\mathrm{d}\ell\mathrm{d} m
 * \f]
 *
 * where \f$ I(\ell, m) \f$ is the model intensity in the sky.
 * For a point source, \f$ I(\ell, m) \f$ is a delta function, with given
 * intensity $\mathrm{I}$, for a Gaussian source oriented along the coordinate
 * axes $\ell$ and $m$, we have
 *
 * \f[
 * I(\ell, m) = \mathrm{I} \frac{1}{\sqrt{2\pi\sigma_\ell\sigma_m}
 *              e^{-\frac{\ell^2}{2\sigma_\ell^2}-\frac{m^2}{2\sigma_m^2}},
 * \]
 *
 * where \f$ \sigma_\ell \f$ and \f$ \sigma_m \f$ are computed from the FWHM
 * major and minor axis in the sky model. The normalization is such that
 * \mathrm{I} represents the integrated flux of the Gaussian source.
 * For Gaussian sources, a position angle or orientation can be
 * given, representing the orientation of the source. Depending on the value
 * of 'OrientationIsAbsolute' in the model component, this is the orientation
 * with respect to the declination axis (if OrientationIsAbsolute is true) or
 * with respect to the m axis for this observation (if OrientationIsAbsolute
 * is false).
 * To compute Gausian sources, the coordinate system is rotated to a new
 * coordinate system where it is oriented along the axes. For the case where
 * 'OrientationIsAbsolute' is true, the $u,v,w$ coordinates are phase-shifted
 * to the position of the source.
 */

class Simulator : public ModelComponentVisitor {
 public:
  /**
   * @brief Construct a new Simulator object
   *
   * @param reference Reference direction (phase center)
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
    DuoMatrix() = default;

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
    size_t itsNRows{0};
  };

  void simulate(const std::shared_ptr<const ModelComponent>& component);

 private:
  void visit(const PointSource& component) override;
  void visit(const GaussianSource& component) override;

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
