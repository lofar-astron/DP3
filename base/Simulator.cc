// Simulator.cc: Compute visibilities for different model components types
// (implementation of ModelComponentVisitor).
//
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later
//
// $Id$

#include <casacore/casa/Arrays/MatrixMath.h>

#include <cmath>
#include <cstddef>

#include <aocommon/constants.h>

#include "Simulator.h"
#include "GaussianSource.h"
#include "PointSource.h"
#include "../steps/PhaseShift.h"

#include "../common/StreamUtil.h"

namespace dp3 {
namespace base {

namespace {

/**
 * Compute station phase shifts.
 *
 * \f[ \mathrm{stationphases}(p) = \frac{2\pi}{c}((u_p, v_p, w_p) \cdot (\ell,
 * m, n)) \f]
 *
 * \f[ \mathrm{phases}(p) = e^{\mathrm{stationphases}(p)} \f]
 *
 * @param nStation Number of stations
 * @param nChannel Number of channels
 * @param lmn LMN coordinates of source, should be length 3
 * @param uvw Station UVW coordinates, matrix of shape (3, nSt)
 * @param freq Channel frequencies, should be length nChannel
 * @param shift Output matrix (2 for real,imag), shift per station, matrix of
 * shape (3, nSt)
 * @param stationPhases Output vector, store per station \f$(x_1,y_1)\f$
 */
void phases(size_t nStation, size_t nChannel, const double* lmn,
            const xt::xtensor<double, 2>& uvw, const std::vector<double>& freq,
            Simulator::DuoMatrix<double>& shift,
            std::vector<double>& scaled_ncp_uvw,
            std::vector<double>& stationPhases,
            std::vector<double>& stationEarthRotation);

float sinc(float x);

void spectrum(const PointSource& component, size_t nChannel,
              const std::vector<double>& freq,
              Simulator::DuoMatrix<double>& spectrum, bool stokesIOnly);
}  // Unnamed namespace.

Simulator::Simulator(const Direction& reference, size_t nStation,
                     const std::vector<Baseline>& baselines,
                     const std::vector<double>& freq,
                     const std::vector<double>& chanWidths,
                     const xt::xtensor<double, 2>& stationUVW,
                     casacore::Cube<dcomplex>& buffer, bool correctFreqSmearing,
                     bool stokesIOnly)
    : itsReference(reference),
      itsNStation(nStation),
      itsNBaseline(baselines.size()),
      itsNChannel(freq.size()),
      itsCorrectTimeSmearing(false),
      itsCorrectFreqSmearing(correctFreqSmearing),
      itsStokesIOnly(stokesIOnly),
      itsBaselines(baselines),
      itsFreq(freq),
      itsChanWidths(chanWidths),
      itsScaledNcpUvw(3, 0.0),
      itsStationUVW(&stationUVW),
      itsBuffer(buffer),
      itsShiftBuffer(),
      itsSpectrumBuffer() {
  itsShiftBuffer.resize(itsNChannel, nStation);
  itsStationPhases.resize(nStation);
  itsStationEarthRotation.resize(nStation);
  if (stokesIOnly) {
    itsSpectrumBuffer.resize(1, itsNChannel);
  } else {
    itsSpectrumBuffer.resize(4, itsNChannel);
  }
}

Simulator::Simulator(const Direction& reference, size_t nStation,
                     const std::vector<Baseline>& baselines,
                     const std::vector<double>& freq,
                     const std::vector<double>& chanWidths,
                     const std::vector<double>& scaled_ncp_uvw,
                     const xt::xtensor<double, 2>& stationUVW,
                     casacore::Cube<dcomplex>& buffer, bool correctTimeSmearing,
                     bool correctFreqSmearing, bool stokesIOnly)
    : Simulator(reference, nStation, baselines, freq, chanWidths, stationUVW,
                buffer, correctFreqSmearing, stokesIOnly) {
  itsScaledNcpUvw = scaled_ncp_uvw;
  itsCorrectTimeSmearing = correctTimeSmearing;
}

void Simulator::simulate(
    const std::shared_ptr<const ModelComponent>& component) {
  component->accept(*this);
}

void Simulator::visit(const PointSource& component) {
  // Compute LMN coordinates.
  double lmn[3];
  radec2lmn(itsReference, component.direction(), lmn);
  // Compute station phase shifts.
  phases(itsNStation, itsNChannel, lmn, *itsStationUVW, itsFreq, itsShiftBuffer,
         itsScaledNcpUvw, itsStationPhases, itsStationEarthRotation);

  // Compute component spectrum.
  spectrum(component, itsNChannel, itsFreq, itsSpectrumBuffer, itsStokesIOnly);

  // Set number of correlations
  int nCorr = 4;
  if (itsStokesIOnly) {
    nCorr = 1;
  }

  // Use preallocated storage (outside for loops)
  std::vector<double> temp_prod_real_I(itsNChannel);
  std::vector<double> temp_prod_imag_I(itsNChannel);
  std::vector<double> smear_terms(itsNChannel);

  std::vector<double> temp_prod_real_Q;
  std::vector<double> temp_prod_imag_Q;
  std::vector<double> temp_prod_real_U;
  std::vector<double> temp_prod_imag_U;
  std::vector<double> temp_prod_real_V;
  std::vector<double> temp_prod_imag_V;

  // If needed, preallocate storage for Q,U,V
  if (!itsStokesIOnly) {
    temp_prod_real_Q.resize(itsNChannel);
    temp_prod_imag_Q.resize(itsNChannel);
    temp_prod_real_U.resize(itsNChannel);
    temp_prod_imag_U.resize(itsNChannel);
    temp_prod_real_V.resize(itsNChannel);
    temp_prod_imag_V.resize(itsNChannel);
  }

  // Compute visibilities.
  // The following loop can be parallelized, but because there is already
  // a parallelization over sources, this is not necessary.
  // It used to be parallelized with a #pragma omp parallel for, but since
  // the outer loop was also parallelized with omp, this had no effect
  // (since omp by default doesn't parallelize nested loops). After the change
  // to a ThreadPool, this would parallize each separate source, which started
  // hundreds of threads on many-core machines.
  for (size_t bl = 0; bl < itsNBaseline; ++bl) {
    dcomplex* buffer = &itsBuffer(0, 0, bl);
    const size_t p = itsBaselines[bl].first;
    const size_t q = itsBaselines[bl].second;

    if (p == q) {
      buffer += itsNChannel * nCorr;
    } else {
      // Note the notation:
      // Each complex number is represented as  (x+ j y)
      // where x: real part, y: imaginary part
      // Subscripts _p and _q denote stations p and q
      // Subscript _c denotes channel (SpectrumBuffer)
      const double* x_p = &(itsShiftBuffer.real(0, p));
      const double* y_p = &(itsShiftBuffer.imag(0, p));
      const double* x_q = &(itsShiftBuffer.real(0, q));
      const double* y_q = &(itsShiftBuffer.imag(0, q));
      const double* x_c = itsSpectrumBuffer.realdata();
      const double* y_c = itsSpectrumBuffer.imagdata();

      // Precompute smearing factors if needed
      if (itsCorrectTimeSmearing) {
#pragma GCC ivdep
        for (size_t ch = 0; ch < itsNChannel; ++ch) {
          // Smearing is the sinc of half the phase change over the integration
          // interval
          smear_terms[ch] = sinc(
              0.5 * (itsStationEarthRotation[q] - itsStationEarthRotation[p]) *
              itsFreq[ch]);
        }
      } else {
        std::fill(smear_terms.begin(), smear_terms.end(), 1.0);
      }
      if (itsCorrectFreqSmearing) {
#pragma GCC ivdep
        for (size_t ch = 0; ch < itsNChannel; ++ch) {
          // Smearing is the sinc of half the phase change over the integration
          // interval
          smear_terms[ch] *=
              sinc(0.5 * (itsStationPhases[q] - itsStationPhases[p]) *
                   itsChanWidths[ch]);
        }
      }

      if (itsStokesIOnly) {
#pragma GCC ivdep
        for (size_t ch = 0; ch < itsNChannel; ++ch) {
          // Compute baseline phase shift.
          // Compute visibilities.
          const double q_conj_p_real = (*x_p) * (*x_q) + (*y_p) * (*y_q);
          const double q_conj_p_imag = (*x_p) * (*y_q) - (*x_q) * (*y_p);
          temp_prod_real_I[ch] =
              q_conj_p_real * (*x_c) - q_conj_p_imag * (*y_c);
          temp_prod_imag_I[ch] =
              q_conj_p_real * (*y_c++) + q_conj_p_imag * (*x_c++);
          ++x_p;
          ++y_p;
          ++x_q;
          ++y_q;
        }  // Channels.
        if (itsCorrectTimeSmearing || itsCorrectFreqSmearing) {
#pragma GCC ivdep
          for (size_t ch = 0; ch < itsNChannel; ++ch) {
            *buffer++ += dcomplex(smear_terms[ch] * temp_prod_real_I[ch],
                                  smear_terms[ch] * temp_prod_imag_I[ch]);
          }
        } else {
#pragma GCC ivdep
          for (size_t ch = 0; ch < itsNChannel; ++ch) {
            *buffer++ += dcomplex(temp_prod_real_I[ch], temp_prod_imag_I[ch]);
          }
        }
      } else {
#pragma GCC ivdep
        for (size_t ch = 0; ch < itsNChannel; ++ch) {
          // Compute baseline phase shift.
          // Compute visibilities.
          const double q_conj_p_real = (*x_p) * (*x_q) + (*y_p) * (*y_q);
          const double q_conj_p_imag = (*x_p) * (*y_q) - (*x_q) * (*y_p);
          temp_prod_real_I[ch] =
              q_conj_p_real * (*x_c) - q_conj_p_imag * (*y_c);
          temp_prod_imag_I[ch] =
              q_conj_p_real * (*y_c++) + q_conj_p_imag * (*x_c++);
          temp_prod_real_Q[ch] =
              q_conj_p_real * (*x_c) - q_conj_p_imag * (*y_c);
          temp_prod_imag_Q[ch] =
              q_conj_p_real * (*y_c++) + q_conj_p_imag * (*x_c++);
          temp_prod_real_U[ch] =
              q_conj_p_real * (*x_c) - q_conj_p_imag * (*y_c);
          temp_prod_imag_U[ch] =
              q_conj_p_real * (*y_c++) + q_conj_p_imag * (*x_c++);
          temp_prod_real_V[ch] =
              q_conj_p_real * (*x_c) - q_conj_p_imag * (*y_c);
          temp_prod_imag_V[ch] =
              q_conj_p_real * (*y_c++) + q_conj_p_imag * (*x_c++);
          ++x_p;
          ++y_p;
          ++x_q;
          ++y_q;
        }  // Channels.

        if (itsCorrectFreqSmearing) {
#pragma GCC ivdep
          for (size_t ch = 0; ch < itsNChannel; ++ch) {
            *buffer++ += dcomplex(smear_terms[ch] * temp_prod_real_I[ch],
                                  smear_terms[ch] * temp_prod_imag_I[ch]);
            *buffer++ += dcomplex(smear_terms[ch] * temp_prod_real_Q[ch],
                                  smear_terms[ch] * temp_prod_imag_Q[ch]);
            *buffer++ += dcomplex(smear_terms[ch] * temp_prod_real_U[ch],
                                  smear_terms[ch] * temp_prod_imag_U[ch]);
            *buffer++ += dcomplex(smear_terms[ch] * temp_prod_real_V[ch],
                                  smear_terms[ch] * temp_prod_imag_V[ch]);
          }
        } else {
#pragma GCC ivdep
          for (size_t ch = 0; ch < itsNChannel; ++ch) {
            *buffer++ += dcomplex(temp_prod_real_I[ch], temp_prod_imag_I[ch]);
            *buffer++ += dcomplex(temp_prod_real_Q[ch], temp_prod_imag_Q[ch]);
            *buffer++ += dcomplex(temp_prod_real_U[ch], temp_prod_imag_U[ch]);
            *buffer++ += dcomplex(temp_prod_real_V[ch], temp_prod_imag_V[ch]);
          }
        }
      }
    }
  }  // Baselines.
}

void Simulator::visit(const GaussianSource& component) {
  // Compute LMN coordinates.
  double lmn[3];
  radec2lmn(itsReference, component.direction(), lmn);

  // Compute station phase shifts.
  phases(itsNStation, itsNChannel, lmn, *itsStationUVW, itsFreq, itsShiftBuffer,
         itsScaledNcpUvw, itsStationPhases, itsStationEarthRotation);

  // Compute component spectrum.
  spectrum(component, itsNChannel, itsFreq, itsSpectrumBuffer, itsStokesIOnly);

  // Convert position angle from North over East to the angle used to
  // rotate the right-handed UV-plane.
  const double phi = M_PI_2 + component.getPositionAngle() + M_PI;
  const double cosPhi = cos(phi);
  const double sinPhi = sin(phi);

  casacore::Matrix<double> uvwShifted;

  // Create a casacore view on itsStationUVW for now.
  const casacore::IPosition stationUvwShape(2, itsStationUVW->shape(1),
                                            itsStationUVW->shape(0));
  // Creating the view unfortunately requires a const_cast.
  const casacore::Matrix<double> casaStationUvw(
      stationUvwShape, const_cast<double*>(itsStationUVW->data()),
      casacore::SHARE);

  if (component.getPositionAngleIsAbsolute()) {
    // Correct for projection and rotation effects: phase shift u, v, w to the
    // position of the source for evaluating the gaussian
    casacore::Matrix<double> euler_matrix_phasecenter(3, 3);
    casacore::Matrix<double> euler_matrix_source(3, 3);
    dp3::steps::PhaseShift::fillEulerMatrix(euler_matrix_phasecenter,
                                            itsReference);
    dp3::steps::PhaseShift::fillEulerMatrix(euler_matrix_source,
                                            component.direction());
    casacore::Matrix<double> euler_matrix = casacore::product(
        casacore::transpose(euler_matrix_source), euler_matrix_phasecenter);

    uvwShifted = casacore::product(euler_matrix, casaStationUvw);
  } else {
    uvwShifted = casaStationUvw;
  }

  // Take care of the conversion of axis lengths from FWHM in radians to
  // sigma.
  const double fwhm2sigma = 1.0 / (2.0 * std::sqrt(2.0 * std::log(2.0)));
  const double uScale = component.getMajorAxis() * fwhm2sigma;
  const double vScale = component.getMinorAxis() * fwhm2sigma;

  // Set number of correlations
  int nCorr = 4;
  if (itsStokesIOnly) {
    nCorr = 1;
  }

  std::vector<double> smear_terms(itsNChannel);
  const double inv_c_sqr = 1.0 / (casacore::C::c * casacore::C::c);

  for (size_t bl = 0; bl < itsNBaseline; ++bl) {
    dcomplex* buffer = &itsBuffer(0, 0, bl);
    const size_t p = itsBaselines[bl].first;
    const size_t q = itsBaselines[bl].second;

    if (p == q) {
      buffer += itsNChannel * nCorr;
    } else {
      double u = uvwShifted(0, q);
      double v = uvwShifted(1, q);

      u -= uvwShifted(0, p);
      v -= uvwShifted(1, p);

      // Rotate (u, v) by the position angle and scale with the major
      // and minor axis lengths (FWHM in rad).
      const double uPrime = uScale * (u * cosPhi - v * sinPhi);
      const double vPrime = vScale * (u * sinPhi + v * cosPhi);

      // Compute uPrime^2 + vPrime^2 and pre-multiply with -2.0 * PI^2
      // / C^2.
      const double uvPrime =
          (-2.0 * M_PI * M_PI) * (uPrime * uPrime + vPrime * vPrime);
      // Note the notation:
      // Each complex number is represented as  (x+ j y)
      // where x: real part, y: imaginary part
      // Subscripts _p and _q denote stations p and q
      // Subscript _c denotes channel (SpectrumBuffer)
      const double* x_p = &(itsShiftBuffer.real(0, p));
      const double* y_p = &(itsShiftBuffer.imag(0, p));
      const double* x_q = &(itsShiftBuffer.real(0, q));
      const double* y_q = &(itsShiftBuffer.imag(0, q));
      const double* x_c = itsSpectrumBuffer.realdata();
      const double* y_c = itsSpectrumBuffer.imagdata();

      // Precompute amplitudes
#pragma GCC ivdep
      for (size_t ch = 0; ch < itsNChannel; ++ch) {
        const double lambda2 = itsFreq[ch] * itsFreq[ch] * inv_c_sqr * uvPrime;
        smear_terms[ch] = exp(lambda2);
      }
      // Precompute smearing factors if needed and modify amplitudes
      if (itsCorrectTimeSmearing) {
#pragma GCC ivdep
        for (size_t ch = 0; ch < itsNChannel; ++ch) {
          // Smearing is the sinc of half the phase change over the integration
          // interval
          smear_terms[ch] *= sinc(
              0.5 * (itsStationEarthRotation[q] - itsStationEarthRotation[p]) *
              itsFreq[ch]);
        }
      }
      if (itsCorrectFreqSmearing) {
#pragma GCC ivdep
        for (size_t ch = 0; ch < itsNChannel; ++ch) {
          // Smearing is the sinc of half the phase change over the integration
          // interval
          smear_terms[ch] *=
              sinc(0.5 * (itsStationPhases[q] - itsStationPhases[p]) *
                   itsChanWidths[ch]);
        }
      }

      if (itsStokesIOnly) {
#pragma GCC ivdep
        for (size_t ch = 0; ch < itsNChannel; ++ch) {
          // Compute baseline phase shift.
          // Compute visibilities.
          const double q_conj_p_real = (*x_p) * (*x_q) + (*y_p) * (*y_q);
          const double q_conj_p_imag = (*x_p) * (*y_q) - (*x_q) * (*y_p);
          const double temp_prod_real_I =
              q_conj_p_real * (*x_c) - q_conj_p_imag * (*y_c);
          const double temp_prod_imag_I =
              q_conj_p_real * (*y_c++) + q_conj_p_imag * (*x_c++);
          ++x_p;
          ++y_p;
          ++x_q;
          ++y_q;
          *buffer++ += dcomplex(smear_terms[ch] * temp_prod_real_I,
                                smear_terms[ch] * temp_prod_imag_I);
        }
      } else {
#pragma GCC ivdep
        for (size_t ch = 0; ch < itsNChannel; ++ch) {
          // Compute baseline phase shift.
          // Compute visibilities.
          const double q_conj_p_real = (*x_p) * (*x_q) + (*y_p) * (*y_q);
          const double q_conj_p_imag = (*x_p) * (*y_q) - (*x_q) * (*y_p);
          const double temp_prod_real_I =
              q_conj_p_real * (*x_c) - q_conj_p_imag * (*y_c);
          const double temp_prod_imag_I =
              q_conj_p_real * (*y_c++) + q_conj_p_imag * (*x_c++);
          const double temp_prod_real_Q =
              q_conj_p_real * (*x_c) - q_conj_p_imag * (*y_c);
          const double temp_prod_imag_Q =
              q_conj_p_real * (*y_c++) + q_conj_p_imag * (*x_c++);
          const double temp_prod_real_U =
              q_conj_p_real * (*x_c) - q_conj_p_imag * (*y_c);
          const double temp_prod_imag_U =
              q_conj_p_real * (*y_c++) + q_conj_p_imag * (*x_c++);
          const double temp_prod_real_V =
              q_conj_p_real * (*x_c) - q_conj_p_imag * (*y_c);
          const double temp_prod_imag_V =
              q_conj_p_real * (*y_c++) + q_conj_p_imag * (*x_c++);
          ++x_p;
          ++y_p;
          ++x_q;
          ++y_q;

          // The following operations are memory-bound, unlike
          // the compute-bound segment above. By merging these together,
          // an improvement in speed is expected
          *buffer++ += dcomplex(smear_terms[ch] * temp_prod_real_I,
                                smear_terms[ch] * temp_prod_imag_I);
          *buffer++ += dcomplex(smear_terms[ch] * temp_prod_real_Q,
                                smear_terms[ch] * temp_prod_imag_Q);
          *buffer++ += dcomplex(smear_terms[ch] * temp_prod_real_U,
                                smear_terms[ch] * temp_prod_imag_U);
          *buffer++ += dcomplex(smear_terms[ch] * temp_prod_real_V,
                                smear_terms[ch] * temp_prod_imag_V);
        }
      }
    }
  }  // Baselines.
}

namespace {

inline float sinc(float x) {
  return ((x == 0.0f) ? 1.0f : std::fabs(std::sin(x) / x));
}

// Compute station phase shifts.
inline void phases(size_t nStation, size_t nChannel, const double* lmn,
                   const xt::xtensor<double, 2>& uvw,
                   const std::vector<double>& freq,
                   Simulator::DuoMatrix<double>& shift,
                   std::vector<double>& scaled_ncp_uvw,
                   std::vector<double>& stationPhases,
                   std::vector<double>& stationEarthRotation) {
  double* shiftdata_re = shift.realdata();
  double* shiftdata_im = shift.imagdata();
  constexpr double cinv = 2.0 * M_PI / aocommon::kSpeedOfLight;
#pragma GCC ivdep
  for (size_t st = 0; st < nStation; ++st) {
    stationPhases[st] = cinv * (uvw(st, 0) * lmn[0] + uvw(st, 1) * lmn[1] +
                                uvw(st, 2) * (lmn[2] - 1.0));
    // See https://wsclean.readthedocs.io/en/latest/time_frequency_smearing.html
    // for an explanation of the equation implemented below
    // Note that scaled_ncp_uvw[0] is currently always zero, so these terms have
    // been commented out. When the initialization of scaled_ncp_uvw[0] is
    // changed to some non-zero value, these terms need to be uncommented
    stationEarthRotation[st] =
        cinv * (lmn[0] * (uvw(st, 1) * scaled_ncp_uvw[2] -
                          uvw(st, 2) * scaled_ncp_uvw[1]) +
                lmn[1] * (/* uvw(st,2) * scaled_ncp_uvw[0] */ -uvw(st, 0) *
                          scaled_ncp_uvw[2]) +
                (lmn[2] - 1.0) *
                    (uvw(st, 0) *
                     scaled_ncp_uvw[1] /* - uvw(st,1) * scaled_ncp_uvw[0] */));
  }

  std::vector<double> phase_terms(nChannel);
  // Note that sincos() does not vectorize yet, and
  // separate sin() cos() is merged to sincos() by the compiler.
  // Hence split the loop into separate sin(), cos() loops.
  for (size_t st = 0; st < nStation; ++st) {
#pragma GCC ivdep
    for (size_t ch = 0; ch < nChannel; ++ch) {
      phase_terms[ch] = stationPhases[st] * freq[ch];
    }  // Channels.
#pragma GCC ivdep
    for (size_t ch = 0; ch < nChannel; ++ch) {
      *shiftdata_im++ = std::sin(phase_terms[ch]);
    }  // Channels.
#pragma GCC ivdep
    for (size_t ch = 0; ch < nChannel; ++ch) {
      *shiftdata_re++ = std::cos(phase_terms[ch]);
    }  // Channels.
  }    // Stations.
}

// Compute component spectrum.
inline void spectrum(const PointSource& component, size_t nChannel,
                     const std::vector<double>& freq,
                     Simulator::DuoMatrix<double>& spectrum,
                     bool stokesIOnly = false) {
#pragma GCC ivdep
  for (size_t ch = 0; ch < nChannel; ++ch) {
    Stokes stokes = component.stokes(freq[ch]);

    if (stokesIOnly) {
      spectrum.real(0, ch) = stokes.I;
      spectrum.imag(0, ch) = 0.0;
    } else {
      spectrum.real(0, ch) = stokes.I + stokes.Q;
      spectrum.imag(0, ch) = 0.0;
      spectrum.real(1, ch) = stokes.U;
      spectrum.imag(1, ch) = stokes.V;
      spectrum.real(2, ch) = stokes.U;
      spectrum.imag(2, ch) = -stokes.V;
      spectrum.real(3, ch) = stokes.I - stokes.Q;
      spectrum.imag(3, ch) = 0.0;
    }
  }
}

}  // Unnamed namespace.

}  // namespace base
}  // namespace dp3
