//# DPInfo.cc: General info about DPPP data processing attributes like averaging
//# Copyright (C) 2010
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
//#
//# @author Ger van Diepen

#include <lofar_config.h>
#include <DPPP/DPInfo.h>
#include <Common/LofarLogger.h>
#include <measures/Measures/MeasConvert.h>
#include <measures/Measures/MCPosition.h>
#include <casa/Arrays/ArrayMath.h>

using namespace casa;

namespace LOFAR {
  namespace DPPP {

    DPInfo::DPInfo()
      : itsNeedVisData  (false),
        itsNeedWrite    (false),
        itsNCorr        (0),
        itsNChan        (0),
        itsChanAvg      (1),
        itsNTime        (0),
        itsTimeAvg      (1),
        itsStartTime    (0),
        itsTimeInterval (0),
        itsPhaseCenterIsOriginal (true)
    {}

    void DPInfo::init (uint ncorr, uint nchan,
                       uint ntime, double startTime, double timeInterval,
                       const string& msName)
    {
      itsNCorr        = ncorr;
      itsNChan        = nchan;
      itsOrigNChan    = nchan;
      itsNTime        = ntime;
      itsStartTime    = startTime;
      itsTimeInterval = timeInterval;
      itsMSName       = msName;
    }

    void DPInfo::set (const Vector<double>& chanFreqs,
                      const Vector<double>& chanWidths,
                      const Vector<double>& resolutions,
                      const Vector<double>& effectiveBW,
                      double totalBW, double refFreq)
    {
      itsChanFreqs.reference  (chanFreqs);
      itsChanWidths.reference (chanWidths);
      if (resolutions.size() == 0) {
        itsResolutions.reference (chanWidths);
      } else {
        itsResolutions.reference (resolutions);
      }
      if (effectiveBW.size() == 0) {
        itsEffectiveBW.reference (chanWidths);
      } else {
        itsEffectiveBW.reference (effectiveBW);
      }
      if (totalBW == 0) {
        itsTotalBW = sum(itsEffectiveBW);
      } else {
        itsTotalBW = totalBW;
      }
      if (refFreq == 0) {
        int n = itsChanFreqs.size();
        // Takes mean of middle elements if n is even; takes middle if odd.
        itsRefFreq = 0.5 * (itsChanFreqs[(n-1)/2] + itsChanFreqs[n/2]);
      } else {
        itsRefFreq = refFreq;
      }
    }

    void DPInfo::set (const MPosition& arrayPos,
                      const MDirection& phaseCenter,
                      const MDirection& delayCenter,
                      const MDirection& tileBeamDir)
    {
      itsArrayPos    = arrayPos;
      itsPhaseCenter = phaseCenter;
      itsDelayCenter = delayCenter;
      itsTileBeamDir = tileBeamDir;
    }

    void DPInfo::set (const Vector<String>& antNames,
                      const vector<MPosition>& antPos,
                      const Vector<Int>& ant1,
                      const Vector<Int>& ant2)
    {
      itsAntNames.reference (antNames);
      itsAntPos = antPos;
      itsAnt1.reference (ant1);
      itsAnt2.reference (ant2);
    }

    void DPInfo::set (const Vector<Int>& ant1,
                      const Vector<Int>& ant2)
    {
      itsAnt1.reference (ant1);
      itsAnt2.reference (ant2);
    }

    uint DPInfo::update (uint chanAvg, uint timeAvg)
    {
      if (chanAvg > itsNChan) {
        chanAvg = itsNChan;
      }
      if (timeAvg > itsNTime) {
        timeAvg = itsNTime;
      }
      ASSERTSTR (itsNChan % chanAvg == 0,
                 "When averaging, nr of channels must divide integrally; "
                 "itsNChan=" << itsNChan << " chanAvg=" << chanAvg);
      itsChanAvg *= chanAvg;
      itsNChan = (itsNChan + chanAvg - 1) / chanAvg;
      itsTimeAvg *= timeAvg;
      itsNTime = (itsNTime + timeAvg - 1) / timeAvg;
      itsTimeInterval *= timeAvg;
      Vector<double> freqs(itsNChan);
      Vector<double> widths(itsNChan, 0.);
      Vector<double> resols(itsNChan, 0.);
      Vector<double> effBWs(itsNChan, 0.);
      double totBW = 0;
      for (uint i=0; i<itsNChan; ++i) {
        freqs[i] = 0.5 * (itsChanFreqs[i*chanAvg] +
                          itsChanFreqs[(i+1)*chanAvg - 1]);
        for (uint j=0; j<chanAvg; ++j) {
          widths[i] += itsChanWidths[i*chanAvg+j];
          resols[i] += itsResolutions[i*chanAvg+j];
          effBWs[i] += itsEffectiveBW[i*chanAvg+j];
        }
        totBW += effBWs[i];
      }
      itsChanFreqs.reference   (freqs);
      itsChanWidths.reference  (widths);
      itsResolutions.reference (resols);
      itsEffectiveBW.reference (effBWs);
      itsTotalBW = totBW;
      return chanAvg;
    }

    void DPInfo::update (uint startChan, uint nchan,
                         const vector<uint>& baselines)
    {
      Slice slice(startChan, nchan);
      itsChanFreqs.reference  (itsChanFreqs (slice).copy());
      itsChanWidths.reference (itsChanWidths(slice).copy());
      itsNChan = nchan;
      // Keep only selected baselines.
      if (! baselines.empty()) {
        Vector<Int> ant1 (baselines.size());
        Vector<Int> ant2 (baselines.size());
        for (uint i=0; i<baselines.size(); ++i) {
          ant1 = itsAnt1[i];
          ant2 = itsAnt2[i];
        }
        itsAnt1.reference (ant1);
        itsAnt2.reference (ant2);
        // Clear; they'll be recalculated if needed.
        itsBLength.resize (0);
        itsAutoCorrIndex.resize (0);
      }
    }

    const vector<double>& DPInfo::getBaselineLengths() const
    {
      // Calculate the baseline lengths if not done yet.
      if (itsBLength.empty()) {
        // First get the antenna positions.
        const vector<MPosition>& antPos = antennaPos();
        vector<Vector<double> > antVec;
        antVec.reserve (antPos.size());
        for (vector<MPosition>::const_iterator iter = antPos.begin();
             iter != antPos.end(); ++iter) {
          // Convert to ITRF and keep as x,y,z in m.
          antVec.push_back
           (MPosition::Convert(*iter, MPosition::ITRF)().getValue().getValue());
        }
        // Fill in the length of each baseline.
        vector<double> blength;
        itsBLength.reserve (itsAnt1.size());
        for (uint i=0; i<itsAnt1.size(); ++i) {
          Array<double> diff(antVec[itsAnt2[i]] - antVec[itsAnt1[i]]);
          itsBLength.push_back (sqrt(sum(diff*diff)));
        }
      }
      return itsBLength;
    }

    const vector<int>& DPInfo::getAutoCorrIndex() const
    {
      if (itsAutoCorrIndex.empty()) {
        int nant = 1 + std::max(max(itsAnt1), max(itsAnt2));
        itsAutoCorrIndex.resize (nant);
        std::fill (itsAutoCorrIndex.begin(), itsAutoCorrIndex.end(), -1);
        // Keep the baseline table index for the autocorrelations.
        for (uint i=0; i<itsAnt1.size(); ++i) {
          if (itsAnt1[i] == itsAnt2[i]) {
            itsAutoCorrIndex[itsAnt1[i]] = i;
          }
        }
      }
      return itsAutoCorrIndex;
    }

  } //# end namespace
}
