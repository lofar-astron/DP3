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

#include "DPInfo.h"
#include "DPInput.h"
#include "Exceptions.h"

#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/MCPosition.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayIO.h>
#include <casacore/casa/BasicSL/STLIO.h>

#include <cassert>

using namespace casacore;
using namespace std;

namespace DP3 {
  namespace DPPP {

    DPInfo::DPInfo()
      : itsNeedVisData  (false),
        itsWriteData    (false),
        itsWriteFlags   (false),
        itsWriteWeights (false),
        itsMetaChanged  (false),
        itsNCorr        (0),
        itsStartChan    (0),
        itsNChan        (0),
        itsChanAvg      (1),
        itsNTime        (0),
        itsTimeAvg      (1),
        itsStartTime    (0),
        itsTimeInterval (0),
        itsPhaseCenterIsOriginal (true),
        itsBeamCorrectionMode(NoBeamCorrection),
        itsNThreads     (0)
    {}

    void DPInfo::init (uint ncorr, uint startChan, uint nchan,
                       uint ntime, double startTime, double timeInterval,
                       const string& msName, const string& antennaSet)
    {
      itsNCorr        = ncorr;
      itsStartChan    = startChan;
      itsNChan        = nchan;
      itsOrigNChan    = nchan;
      itsNTime        = ntime;
      itsStartTime    = startTime;
      itsTimeInterval = timeInterval;
      itsMSName       = msName;
      itsAntennaSet   = antennaSet;
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

    void DPInfo::set (const Vector<casacore::String>& antNames,
                      const Vector<Double>& antDiam,
                      const vector<MPosition>& antPos,
                      const Vector<Int>& ant1,
                      const Vector<Int>& ant2)
    {
      assert (antNames.size() == antDiam.size()  &&
              antNames.size() == antPos.size());
      assert (ant1.size() == ant2.size());
      itsAntNames.reference (antNames);
      itsAntDiam.reference (antDiam);
      itsAntPos = antPos;
      itsAnt1.reference (ant1);
      itsAnt2.reference (ant2);
      // Set which antennae are used.
      setAntUsed();
    }

    void DPInfo::setAntUsed()
    {
      itsAntUsed.clear();
      itsAntMap.resize (itsAntNames.size());
      std::fill (itsAntMap.begin(), itsAntMap.end(), -1);
      for (uint i=0; i<itsAnt1.size(); ++i) {
        assert (itsAnt1[i] >= 0  &&  itsAnt1[i] < int(itsAntMap.size())  &&
                itsAnt2[i] >= 0  &&  itsAnt2[i] < int(itsAntMap.size()));
        itsAntMap[itsAnt1[i]] = 0;
        itsAntMap[itsAnt2[i]] = 0;
      }
      itsAntUsed.reserve (itsAntNames.size());
      for (uint i=0; i<itsAntMap.size(); ++i) {
        if (itsAntMap[i] == 0) {
          itsAntMap[i] = itsAntUsed.size();
          itsAntUsed.push_back (i);
        }
      }
    }

    MeasureHolder DPInfo::copyMeasure(const MeasureHolder fromMeas)
    {
      Record rec;
      String msg;
      assert (fromMeas.toRecord (msg, rec));
      MeasureHolder mh2;
      assert (mh2.fromRecord (msg, rec));
      return mh2;
    }

    uint DPInfo::update (uint chanAvg, uint timeAvg)
    {
      if (chanAvg > itsNChan) {
        chanAvg = itsNChan;
      }
      if (timeAvg > itsNTime) {
        timeAvg = itsNTime;
      }
      if (itsNChan % chanAvg != 0)
        throw Exception("When averaging, nr of channels must divide integrally; "
                        "itsNChan=" + std::to_string(itsNChan) + " chanAvg=" + std::to_string(chanAvg));
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
                         const vector<uint>& baselines, bool removeAnt)
    {
      Slice slice(startChan, nchan);
      itsStartChan=startChan;
      itsChanFreqs.reference  (itsChanFreqs (slice).copy());
      itsChanWidths.reference (itsChanWidths(slice).copy());
      itsResolutions.reference (itsResolutions(slice).copy());
      itsEffectiveBW.reference (itsEffectiveBW(slice).copy());
      itsNChan = nchan;
      // Keep only selected baselines.
      if (! baselines.empty()) {
        Vector<Int> ant1 (baselines.size());
        Vector<Int> ant2 (baselines.size());
        for (uint i=0; i<baselines.size(); ++i) {
          ant1[i] = itsAnt1[baselines[i]];
          ant2[i] = itsAnt2[baselines[i]];
        }
        itsAnt1.reference (ant1);
        itsAnt2.reference (ant2);
        // Clear; they'll be recalculated if needed.
        itsBLength.resize (0);
        itsAutoCorrIndex.resize (0);
      }
      setAntUsed();
      // If needed, remove the stations and renumber the baselines.
      if (removeAnt) {
        removeUnusedAnt();
      }
    }

    void DPInfo::removeUnusedAnt()
    {
      if (itsAntUsed.size() < itsAntMap.size()) {
        // First remove stations.
        Vector<casacore::String> antNames (itsAntUsed.size());
        Vector<Double> antDiam (itsAntUsed.size());
        vector<MPosition> antPos;
        antPos.reserve (itsAntUsed.size());
        for (uint i=0; i<itsAntUsed.size(); ++i) {
          antNames[i] = itsAntNames[itsAntUsed[i]];
          antDiam[i]  = itsAntDiam[itsAntUsed[i]];
          antPos.push_back (itsAntPos[itsAntUsed[i]]);
        }
        // Use the new vectors.
        itsAntNames.reference (antNames);
        itsAntDiam.reference (antDiam);
        itsAntPos.swap (antPos);
        // Renumber the baselines.
        for (uint i=0; i<itsAnt1.size(); ++i) {
          itsAnt1[i] = itsAntMap[itsAnt1[i]];
          itsAnt2[i] = itsAntMap[itsAnt2[i]];
        }
        // Now fill the itsAntUsed and itsAntMap vectors again.
        setAntUsed();
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

    Record DPInfo::toRecord() const
    {
      Record rec;
      rec.define ("NeedVisData", itsNeedVisData);
      rec.define ("WriteData", itsWriteData);
      rec.define ("WriteFlags", itsWriteFlags);
      rec.define ("WriteWeights", itsWriteWeights);
      rec.define ("MetaChanged", itsMetaChanged);
      rec.define ("MSName", itsMSName);
      rec.define ("AntennaSet", itsAntennaSet);
      rec.define ("NCorr", itsNCorr);
      rec.define ("StartChan", itsStartChan);
      rec.define ("OrigNChan", itsOrigNChan);
      rec.define ("NChan", itsNChan);
      rec.define ("ChanAvg", itsChanAvg);
      rec.define ("NTime", itsNTime);
      rec.define ("TimeAvg", itsTimeAvg);
      rec.define ("StartTime", itsStartTime);
      rec.define ("TimeInterval", itsTimeInterval);
      rec.define ("ChanFreqs", itsChanFreqs);
      rec.define ("ChanWidths", itsChanWidths);
      rec.define ("Resolutions", itsResolutions);
      rec.define ("EffectiveBW", itsEffectiveBW);
      rec.define ("TotalBW", itsTotalBW);
      rec.define ("RefFreq", itsRefFreq);
      rec.define ("AntNames", itsAntNames);
      rec.define ("AntDiam", itsAntDiam);
      rec.define ("AntUsed", Vector<int>(itsAntUsed));
      rec.define ("AntMap", Vector<int>(itsAntMap));
      rec.define ("Ant1", itsAnt1);
      rec.define ("Ant2", itsAnt2);
      rec.define ("BLength", Vector<double>(itsBLength));
      rec.define ("AutoCorrIndex", Vector<int>(itsAutoCorrIndex));
      return rec;
    }

    void DPInfo::fromRecord (const Record& rec)
    {
      if (rec.isDefined ("NeedVisData")) {
        rec.get ("NeedVisData", itsNeedVisData);
      }
      if (rec.isDefined ("WriteData")) {
        rec.get ("WriteData", itsWriteData);
      }
      if (rec.isDefined ("WriteFlags")) {
        rec.get ("WriteFlags", itsWriteFlags);
      }
      if (rec.isDefined ("WriteWeights")) {
        rec.get ("WriteWeights", itsWriteWeights);
      }
      if (rec.isDefined ("MetaChanged")) {
        rec.get ("MetaChanged", itsMetaChanged);
      }
      if (rec.isDefined ("MSName")) {
        itsMSName = rec.asString ("MSName");
      }
      if (rec.isDefined ("AntennaSet")) {
        itsAntennaSet = rec.asString ("AntennaSet");
      }
      if (rec.isDefined ("NCorr")) {
        rec.get ("NCorr", itsNCorr);
      }
      if (rec.isDefined ("StartChan")) {
        rec.get ("StartChan", itsStartChan);
      }
      if (rec.isDefined ("OrigNChan")) {
        rec.get ("OrigNChan", itsOrigNChan);
      }
      if (rec.isDefined ("NChan")) {
        rec.get ("NChan", itsNChan);
      }
      if (rec.isDefined ("ChanAvg")) {
        rec.get ("ChanAvg", itsChanAvg);
      }
      if (rec.isDefined ("NTime")) {
        rec.get ("NTime", itsNTime);
      }
      if (rec.isDefined ("TimeAvg")) {
        rec.get ("TimeAvg", itsTimeAvg);
      }
      if (rec.isDefined ("StartTime")) {
        rec.get ("StartTime", itsStartTime);
      }
      if (rec.isDefined ("TimeInterval")) {
        rec.get ("TimeInterval", itsTimeInterval);
      }
      if (rec.isDefined ("ChanFreqs")) {
        rec.get ("ChanFreqs", itsChanFreqs);
      }
      if (rec.isDefined ("ChanWidths")) {
        rec.get ("ChanWidths", itsChanWidths);
      }
      if (rec.isDefined ("Resolutions")) {
        rec.get ("Resolutions", itsResolutions);
      }
      if (rec.isDefined ("EffectiveBW")) {
        rec.get ("EffectiveBW", itsEffectiveBW);
      }
      if (rec.isDefined ("TotalBW")) {
        rec.get ("TotalBW", itsTotalBW);
      }
      if (rec.isDefined ("RefFreq")) {
        rec.get ("RefFreq", itsRefFreq);
      }
      if (rec.isDefined ("AntNames")) {
        rec.get ("AntNames", itsAntNames);
      }
      if (rec.isDefined ("AntDiam")) {
        rec.get ("AntDiam", itsAntDiam);
      }
      ///if (rec.isDefined ("AntUsed")) {
      ///itsAntUsed = rec.toArrayInt("AntUsed").tovector();
      ///}
      ///if (rec.isDefined ("AntMap")) {
      ///  itsAntMap = rec.toArrayInt("AntMap").tovector();
      ///}
      if (rec.isDefined ("Ant1")) {
        rec.get ("Ant1", itsAnt1);
      }
      if (rec.isDefined ("Ant2")) {
        rec.get ("Ant2", itsAnt2);
      }
      ///if (rec.isDefined ("BLength")) {
      ///  itsBLength = rec.toArrayDouble("BLength").tovector();
      ///}
      ///if (rec.isDefined ("AutoCorrIndex")) {
      ///  itsAutoCorrIndex = rec.toArrayInt("AutoCorrIndex").tovector();
      ///}
    }

  } //# end namespace
}
