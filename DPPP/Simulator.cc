//# Simulator.cc: Compute visibilities for different model components types
//# (implementation of ModelComponentVisitor).
//#
//# Copyright (C) 2012
//# ASTRON (Netherlands Institute for Radio Astronomy)
//# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//#
//# This file is part of the LOFAR software suite.
//# The LOFAR software suite is free software: you can redistribute it and/or
//# modify it under the terms of the GNU General Public License as published
//# by the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The LOFAR software suite is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//#
//# $Id$

#include "Simulator.h"
#include "GaussianSource.h"
#include "PointSource.h"

#include <casacore/casa/BasicSL/Constants.h>

#include "../Common/StreamUtil.h" ///

namespace DP3
{
namespace DPPP
{

namespace
{
// Compute LMN coordinates of \p position relative to \p reference.
//
// \param[in]   reference
// Reference position on the celestial sphere.
// \param[in]   position
// Position of interest on the celestial sphere.
// \param[in]   lmn
// Pointer to a buffer of (at least) length three into which the computed LMN
// coordinates will be written.
void radec2lmn(const Position &reference, const Position &position,
               double* lmn);

void phases(size_t nStation, size_t nChannel, const double* lmn,
            const casacore::Matrix<double>& uvw, const casacore::Vector<double>& freq,
            Simulator::Matrix<dcomplex>& shift);

void spectrum(const PointSource &component, size_t nChannel,
              const casacore::Vector<double>& freq,
              Simulator::Matrix<dcomplex>& spectrum, bool stokesIOnly);
} // Unnamed namespace.

Simulator::Simulator(const Position &reference, size_t nStation,
    size_t nBaseline, size_t nChannel, const casacore::Vector<Baseline>& baselines,
    const casacore::Vector<double>& freq, const casacore::Matrix<double>& uvw,
    casacore::Cube<dcomplex>& buffer, bool stokesIOnly)
    :   itsReference(reference),
        itsNStation(nStation),
        itsNBaseline(nBaseline),
        itsNChannel(nChannel),
        itsStokesIOnly(stokesIOnly),
        itsBaselines(baselines),
        itsFreq(freq),
        itsUVW(uvw),
        itsBuffer(buffer),
        itsShiftBuffer(),
        itsSpectrumBuffer()
{
  itsShiftBuffer.resize(nChannel,nStation);
  if (stokesIOnly) {
    itsSpectrumBuffer.resize(1,nChannel);
  } else {
    itsSpectrumBuffer.resize(4,nChannel);
  }
}

void Simulator::simulate(const ModelComponent::ConstPtr &component)
{
    component->accept(*this);
}


void Simulator::visit(const PointSource &component) {
  // Compute LMN coordinates.
  double lmn[3];
  radec2lmn(itsReference, component.position(), lmn);

  // Compute station phase shifts.
  phases(itsNStation, itsNChannel, lmn, itsUVW, itsFreq, itsShiftBuffer);

  // Compute component spectrum.
  spectrum(component, itsNChannel, itsFreq, itsSpectrumBuffer, itsStokesIOnly);

  // Set number of correlations
  int nCorr = 4;
  if (itsStokesIOnly) {
    nCorr = 1;
  }

  // Compute visibilities.
  // The following loop can be parallelized, but because there is already
  // a parallelization over sources, this is not necessary.
  // It used to be parallelized with a #pragma omp parallel for, but since
  // the outer loop was also parallelized with omp, this had no effect
  // (since omp by default doesn't parallelize nested loops). After the change
  // to a ThreadPool, this would parallize each separate source, which started
  // hundreds of threads on many-core machines.
  for(size_t bl = 0; bl < itsNBaseline; ++bl) {
    dcomplex* buffer=&itsBuffer(0,0,bl);
    const size_t p = itsBaselines[bl].first;
    const size_t q = itsBaselines[bl].second;

    if(p == q) {
      buffer+=itsNChannel*nCorr;
    } else {
      const dcomplex *shiftP = &(itsShiftBuffer(0,p));
      const dcomplex *shiftQ = &(itsShiftBuffer(0,q));
      const dcomplex *spectrum = itsSpectrumBuffer.data();

      if (itsStokesIOnly) {
        for (size_t ch = 0; ch < itsNChannel; ++ch)
        {
          // Compute baseline phase shift.
          // Compute visibilities.
          *buffer++ += (*shiftQ) * conj(*shiftP) * (*spectrum++);
          ++shiftP;
          ++shiftQ;
        } // Channels.
      } else {
        for(size_t ch = 0; ch < itsNChannel; ++ch)
        {
          // Compute baseline phase shift.
          const dcomplex blShift = (*shiftQ) * conj(*shiftP);
          ++shiftP;
          ++shiftQ;

          // Compute visibilities.
          *buffer++ += blShift * (*spectrum++);
          *buffer++ += blShift * (*spectrum++);
          *buffer++ += blShift * (*spectrum++);
          *buffer++ += blShift * (*spectrum++);
        } // Channels.
      }
    }
  } // Baselines.
}

void Simulator::visit(const GaussianSource &component)
{
    // Compute LMN coordinates.
    double lmn[3];
    radec2lmn(itsReference, component.position(), lmn);

    // Compute station phase shifts.
    phases(itsNStation, itsNChannel, lmn, itsUVW, itsFreq, itsShiftBuffer);

    // Compute component spectrum.
    spectrum(component, itsNChannel, itsFreq, itsSpectrumBuffer, itsStokesIOnly);

    // Convert position angle from North over East to the angle used to
    // rotate the right-handed UV-plane.
    // TODO: Can probably optimize by changing the rotation matrix instead.
    const double phi = casacore::C::pi_2 + component.positionAngle() + casacore::C::pi;
    const double cosPhi = cos(phi);
    const double sinPhi = sin(phi);

    // Take care of the conversion of axis lengths from FWHM in radians to
    // sigma.
    // TODO: Shouldn't the projection from the celestial sphere to the
    // UV-plane be taken into account here?
    const double fwhm2sigma = 1.0 / (2.0 * std::sqrt(2.0 * std::log(2.0)));
    const double uScale = component.majorAxis() * fwhm2sigma;
    const double vScale = component.minorAxis() * fwhm2sigma;

    // Set number of correlations
    int nCorr = 4;
    if (itsStokesIOnly) {
      nCorr = 1;
    }

    for(size_t bl = 0; bl < itsNBaseline; ++bl)
    {
        dcomplex* buffer=&itsBuffer(0,0,bl);
        const size_t p = itsBaselines[bl].first;
        const size_t q = itsBaselines[bl].second;

        if(p == q) {
          buffer+=itsNChannel*nCorr;
        } else {
            double u = itsUVW(0,q);
            double v = itsUVW(1,q);

            u -= itsUVW(0,p);
            v -= itsUVW(1,p);

            // Rotate (u, v) by the position angle and scale with the major
            // and minor axis lengths (FWHM in rad).
            const double uPrime = uScale * (u * cosPhi - v * sinPhi);
            const double vPrime = vScale * (u * sinPhi + v * cosPhi);

            // Compute uPrime^2 + vPrime^2 and pre-multiply with -2.0 * PI^2
            // / C^2.
            const double uvPrime = (-2.0 * casacore::C::pi * casacore::C::pi)
                * (uPrime * uPrime + vPrime * vPrime);

            const dcomplex *shiftP = &(itsShiftBuffer(0,p));
            const dcomplex *shiftQ = &(itsShiftBuffer(0,q));
            const dcomplex *spectrum = itsSpectrumBuffer.data();

            if (itsStokesIOnly) {
                for (size_t ch = 0; ch < itsNChannel; ++ch)
                {
                    // Compute baseline phase shift.
                    dcomplex blShift = (*shiftQ) * conj(*shiftP);
                    ++shiftP;
                    ++shiftQ;

                    const double ampl = exp((itsFreq[ch] * itsFreq[ch])
                        / (casacore::C::c * casacore::C::c) * uvPrime);

                    blShift *= ampl;

                    // Compute visibilities.
                    *buffer++ += blShift * (*spectrum++);
                } // Channels.
            } else {
                for(size_t ch = 0; ch < itsNChannel; ++ch)
                {
                    // Compute baseline phase shift.
                    dcomplex blShift = (*shiftQ) * conj(*shiftP);
                    ++shiftP;
                    ++shiftQ;

                    const double ampl = exp((itsFreq[ch] * itsFreq[ch])
                        / (casacore::C::c * casacore::C::c) * uvPrime);

                    blShift *= ampl;

                    // Compute visibilities.
                    *buffer++ += blShift * (*spectrum++);
                    *buffer++ += blShift * (*spectrum++);
                    *buffer++ += blShift * (*spectrum++);
                    *buffer++ += blShift * (*spectrum++);
                } // Channels.
            }
        }
    } // Baselines.
}

namespace
{
inline void radec2lmn(const Position &reference, const Position &position,
                      double* lmn)
{
  const double dRA = position[0] - reference[0];
  const double pDEC = position[1];
  const double rDEC = reference[1];
  const double cDEC = cos(pDEC);

  const double l = cDEC * sin(dRA);
  const double m = sin(pDEC) * cos(rDEC) - cDEC * sin(rDEC) * cos(dRA);

  lmn[0] = l;
  lmn[1] = m;
  lmn[2] = sqrt(1.0 - l * l - m * m);
}

// Compute station phase shifts.
inline void phases(size_t nStation, size_t nChannel, const double* lmn,
                   const casacore::Matrix<double>& uvw,
                   const casacore::Vector<double>& freq,
                   Simulator::Matrix<dcomplex>& shift)
{
    dcomplex* shiftdata=shift.data();
    for(size_t st = 0; st < nStation; ++st)
    {
        const double phase = casacore::C::_2pi * (uvw(0,st) * lmn[0]
            + uvw(1,st) * lmn[1] + uvw(2,st) * (lmn[2] - 1.0));

        for(size_t ch = 0; ch < nChannel; ++ch)
        {
            const double chPhase = phase * freq[ch] / casacore::C::c;
            *shiftdata = dcomplex(cos(chPhase), sin(chPhase));
            ++shiftdata;
        } // Channels.
    } // Stations.
}


// Compute component spectrum.
inline void spectrum(const PointSource &component, size_t nChannel,
                     const casacore::Vector<double>& freq,
                     Simulator::Matrix<dcomplex>& spectrum, bool stokesIOnly=false)
{
    for(size_t ch = 0; ch < nChannel; ++ch)
    {
        Stokes stokes = component.stokes(freq[ch]);

        if (stokesIOnly) {
          spectrum(0,ch) = dcomplex(stokes.I, 0.0);
        } else {
          spectrum(0,ch) = dcomplex(stokes.I + stokes.Q, 0.0);
          spectrum(1,ch) = dcomplex(stokes.U, stokes.V);
          spectrum(2,ch) = dcomplex(stokes.U, -stokes.V);
          spectrum(3,ch) = dcomplex(stokes.I - stokes.Q, 0.0);
        }
    }
}
} // Unnamed namespace.

} //# namespace DPPP
} //# namespace LOFAR
