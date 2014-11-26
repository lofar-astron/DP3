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

#include <lofar_config.h>
#include <DPPP/Simulator.h>
#include <DPPP/GaussianSource.h>
#include <DPPP/PointSource.h>
#include <casa/BasicSL/Constants.h>
#include <Common/StreamUtil.h> ///

namespace LOFAR
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

void spectrum(const PointSource &component, size_t nChannel,
    vector<double> freq, dcomplex *spectrum);
} // Unnamed namespace.

Simulator::Simulator(const Position &reference, size_t nStation,
    size_t nBaseline, size_t nChannel, const vector<Baseline> baselines,
    const vector<double> freq, const casa::Matrix<double> uvw,
    casa::Matrix<dcomplex> buffer)
    :   itsReference(reference),
        itsNStation(nStation),
        itsNBaseline(nBaseline),
        itsNChannel(nChannel),
        itsBaselines(baselines),
        itsFreq(freq),
        itsUVW(uvw),
        itsBuffer(buffer),
        itsShiftBuffer(nStation * nChannel),
        itsSpectrumBuffer(nChannel * 4)
{
}

void Simulator::simulate(const ModelComponent::ConstPtr &component)
{
    component->accept(*this);
}


void Simulator::visit(const PointSource &component)
{
    // Compute LMN coordinates.
    double lmn[3];
    radec2lmn(itsReference, component.position(), lmn);
    ///    cout<<"pos="<<itsReference[0]<<' '<<itsReference[1]<<' '<<component.position()[0]<<' '<<component.position()[1]<<endl;
    ///    cout<<"lmn="<<lmn[0]<<' '<<lmn[1]<<' '<<lmn[2]<<endl;

    void phases(size_t nStation, size_t nChannel, double* lmn,
        const casa::Matrix<double> uvw, const vector<double> freq,
        dcomplex* shift);

    // Compute station phase shifts.
    phases(itsNStation, itsNChannel, lmn, itsUVW, itsFreq, &itsShiftBuffer[0]);

    // Compute component spectrum.
    spectrum(component, itsNChannel, itsFreq, &itsSpectrumBuffer[0]);

    // Compute visibilities.
    for(size_t bl = 0; bl < itsNBaseline; ++bl)
    {
        const size_t p = itsBaselines[bl].first;
        const size_t q = itsBaselines[bl].second;

        if(p != q)
        {
            const dcomplex *shiftP = &(itsShiftBuffer[p * itsNChannel]);
            const dcomplex *shiftQ = &(itsShiftBuffer[q * itsNChannel]);
            const dcomplex *spectrum = &(itsSpectrumBuffer[0]);
            for(size_t ch = 0; ch < itsNChannel; ++ch)
            {
                // Compute baseline phase shift.
                const dcomplex blShift = (*shiftQ) * conj(*shiftP);
                ++shiftP;
                ++shiftQ;

                // Compute visibilities.
                *itsBuffer += blShift * (*spectrum);
                ++itsBuffer;
                ++spectrum;
                *itsBuffer += blShift * (*spectrum);
                ++itsBuffer;
                ++spectrum;
                *itsBuffer += blShift * (*spectrum);
                ++itsBuffer;
                ++spectrum;
                *itsBuffer += blShift * (*spectrum);
                ++itsBuffer;
                ++spectrum;

                // Move to next channel.
                itsBuffer -= 4;
                itsBuffer.forward(1);
            } // Channels.

            itsBuffer.backward(1, itsNChannel);
        }

        // Move to next baseline.
        itsBuffer.forward(2);
        ++itsBaselines;
    } // Baselines.

    itsBuffer.backward(2, itsNBaseline);
    itsBaselines -= itsNBaseline;
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

inline void phases(size_t nStation, size_t nChannel, double* lmn,
    const casa::Matrix<double> uvw, const vector<double> freq, dcomplex* shift)
{
    // Compute station phase shifts.
    for(size_t st = 0; st < nStation; ++st)
    {
        const double phase = casa::C::_2pi * (
            uvw[nStation,0] * lmn[0] +
            uvw[nStation,1] * lmn[1] +
            uvw[nStation,2] * (lmn[2] - 1.0));

        for(size_t ch = 0; ch < nChannel; ++ch)
        {
            const double chPhase = phase * freq[ch] / casa::C::c;
            shift[nChannel*st+ch] = dcomplex(cos(chPhase), sin(chPhase));
        } // Channels.
    } // Stations.
}

inline void spectrum(const PointSource &component, size_t nChannel,
    double *freq, dcomplex *spectrum)
{
    // Compute component spectrum.
    for(size_t ch = 0; ch < nChannel; ++ch)
    {
        Stokes stokes = component.stokes(freq[ch]);

        spectrum[4*ch+0] = dcomplex(stokes.I + stokes.Q, 0.0);
        spectrum[4*ch+1] = dcomplex(stokes.U, stokes.V);
        spectrum[4*ch+2] = dcomplex(stokes.U, -stokes.V);
        spectrum[4*ch+3] = dcomplex(stokes.I - stokes.Q, 0.0);
    }
}

} // Unnamed namespace.

} //# namespace DPPP
} //# namespace LOFAR
